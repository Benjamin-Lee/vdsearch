import hashlib
import logging
from pathlib import Path
from typing import Optional

import nimporter
from numpy import product
import pandas as pd
import skbio
import typer
import rich.progress
from vdsearch.nim import write_seqs as ws
from vdsearch.rich_wrapper import MyTyper
from vdsearch.types import FASTA
from vdsearch.utils import typer_unpacker

app = MyTyper(hidden=True)


@app.command()
@typer_unpacker
def dbn2tsv(
    dbn: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Dot bracket notation structure file",
    ),
    outfile: Optional[Path] = typer.Option(
        None,
        "--outfile",
        "-o",
        dir_okay=False,
        file_okay=True,
        help="Path to output file",
    ),
):
    """
    Convert an RNAfold .dbn file to a .tsv file.

    If no output file is specified, nothing is written to disk.
    Why? So that the function can by used in a notebook directly.
    """
    result = []

    seq_id, structure, seq, mfe = "", "", "", 0.0
    with dbn.open() as f:
        for i, line in enumerate(f):
            if i % 3 == 0:
                seq_id = line.strip()[1:]
            elif i % 3 == 1:
                seq = line.strip()
            elif i % 3 == 2:
                structure = line.strip()[0 : line.strip().rfind("(")]
                mfe = float(line.strip().split("(")[-1].split(")")[0])
                paired_percent = (len(structure) - structure.count(".")) / len(
                    structure
                )

                # count how many hairpins there are in the structure
                # a hairpin is defined by going from a ( to a )
                hairpins = 0
                last_char = "("
                for char in structure:
                    if char == "(" and last_char == "(":
                        continue
                    elif char == ")" and last_char == ")":
                        continue
                    elif char == ")" and last_char == "(":
                        hairpins += 1
                        last_char = char
                    elif char == "(" and last_char == ")":
                        last_char = char

                # the i + 1 / 3 is to account for the three rows per sequence
                # the - 1 is to make it zero indexed
                result.append(
                    pd.DataFrame(
                        [
                            [
                                seq_id.split()[
                                    0
                                ],  # split on space to get the sequence id
                                structure,
                                seq,
                                mfe,
                                paired_percent,
                                hairpins,
                            ]
                        ],
                        index=[((i + 1) / 3) - 1],
                        columns=[
                            "seq_id",
                            "structure",
                            "seq",
                            "mfe",
                            "paired_percent",
                            "hairpins",
                        ],
                    )
                )

    # either write to disk or return the result, depending on the arguments
    result_df: pd.DataFrame = pd.concat(result)
    if outfile is not None:
        result_df.to_csv(outfile, sep="\t", index=False)
    else:
        return result_df


@app.command()
def cache_path():
    """
    Return the path to the cache directory.
    """
    return Path(typer.get_app_dir("vdsearch"))


@app.command()
def clusters2fastas(
    cluster_tsv: Path = typer.Argument(
        ..., help="Path to cluster TSV file from MMseqs2"
    ),
    fasta: Path = FASTA,
    outdir: Path = typer.Argument(
        ..., dir_okay=True, file_okay=False, help="Path to output directory"
    ),
):
    """
    Write each cluster to a FASTA files.

    This is useful for doing per-cluster analysis with tools that expect a single FASTA file (_e.g._ alignment).
    """

    # map each sequence to its cluster
    seq_to_cluster = {}
    cluster_df = pd.read_csv(cluster_tsv, sep="\t", names=["cluster_id", "seq_id"])
    for cluster_id, seq_id in cluster_df.itertuples(index=False):
        seq_to_cluster[seq_id] = cluster_id

    clusters = cluster_df["cluster_id"].unique().tolist()

    outdir.mkdir(exist_ok=True, parents=True)

    # do a call to Nim to write the FASTA files out
    ws.write_clusters(str(fasta), str(outdir), seq_to_cluster, clusters)

    logging.done(f"Wrote {len(clusters)} cluster FASTA files to {outdir}")  # type: ignore


# rank sequences by their multiplied ribozyme E values
@app.command()
def rank_by_ribozyme(
    infernal_tsv: Path = typer.Argument(..., help="Path to infernal TSV file"),
    outfile: Path = typer.Argument(
        ...,
        dir_okay=False,
        file_okay=True,
        help="Path to output file",
    ),
):
    """
    Rank sequences by their multiplied ribozyme E values.
    """
    df = pd.read_csv(infernal_tsv, sep="\t")
    outdf = []
    for seq_name, seq_df in df.groupby(["seq_id"]):

        e_values = []

        for pol in ["+", "-"]:
            e_values.extend(
                seq_df.query(f"strand == '{pol}'")
                .sort_values(by=["evalue"], ascending=True)
                .drop_duplicates(subset=["seq_id"], keep="first")
                .evalue.values
            )

        outdf.append(dict(seq_id=seq_name, evalue=product(e_values)))
    pd.DataFrame(outdf, columns=["seq_id", "evalue"]).sort_values(
        by=["evalue"], ascending=True
    ).to_csv(outfile, sep="\t", index=False)


@app.command()
def summarize(
    fasta: Path = FASTA,
    ribozyme_tsv: Path = typer.Argument(..., help="Path to ribozyme TSV file"),
    dbn: Path = typer.Argument(..., help="Path to structure dbn file"),
    # raw_fasta: Path = FASTA,
    source: str = typer.Argument(..., help="Source of the sequences"),
    outfile: Path = typer.Argument(..., help="Path to output file"),
    header: bool = typer.Option(
        False, help="Include a header. Only useful for the first run"
    ),
):
    """
    Summarize the results of the analysis.
    """

    sequences = {}
    seq_id_to_new_id = {}
    ribozymes_df = pd.read_csv(ribozyme_tsv, sep="\t")
    dbn_df = dbn2tsv(dbn)

    seq: skbio.DNA
    with rich.progress.open(fasta) as f:
        for seq in skbio.read(f, format="fasta", constructor=skbio.DNA):
            seq_data = {}

            seq_id = seq.metadata["id"]
            # basic sequence information
            seq_data["seq_id"] = seq_id
            seq_data["description"] = seq.metadata["description"]
            seq_data["unit_length"] = len(seq)
            seq_data["gc_content"] = seq.gc_content()

            seq_ribozymes = ribozymes_df.query(f"seq_id == '{seq_id}'").sort_values(
                by=["evalue"], ascending=True
            )
            seq_data["has_ribozymes"] = (
                bool(len(seq_ribozymes))
                and ["Pospi_RY"] != seq_ribozymes.ribozyme.unique().tolist()
            )

            # I've renamed Polarity to symmetry (and made it boolean) but I'm keeping backwards compatibility
            if "Polarity" in seq_ribozymes.columns:
                seq_data["symmetric"] = (
                    seq_ribozymes.Polarity.unique()[0] == "(+) and (-)"
                )
            else:
                seq_data["symmetric"] = seq_ribozymes.symmetric.unique()[0]

            rzs_present = [
                seq_ribozymes.loc[seq_ribozymes.strand == strand]
                for strand in ["+", "-"]
            ]
            # print(rzs_present)
            if rzs_present[0].shape[0]:
                seq_data["rz_plus"] = rzs_present[0].ribozyme.values[0]
                seq_data["rz_plus_evalue"] = rzs_present[0].evalue.values[0]
                seq_data["rz_plus_from"] = rzs_present[0]["from"].values[0]
                seq_data["rz_plus_to"] = rzs_present[0]["to"].values[0]
            if rzs_present[1].shape[0]:
                seq_data["rz_minus"] = rzs_present[1].ribozyme.values[0]
                seq_data["rz_minus_evalue"] = rzs_present[1].evalue.values[0]
                seq_data["rz_minus_from"] = rzs_present[1]["from"].values[0]
                seq_data["rz_minus_to"] = rzs_present[1]["to"].values[0]

            # merge in
            seq_data.update(
                dbn_df.query(f"seq_id == '{seq_id}'").to_dict(orient="records")[0]
            )

            seq_data["source"] = source

            new_id = (
                "NV_"
                + hashlib.blake2b(str(seq).encode("utf-8"), digest_size=8).hexdigest()
            )
            seq_data["vdsearch_id"] = new_id
            seq_id_to_new_id[seq_id] = new_id
            sequences[new_id] = seq_data

    # find the original lengths
    # with Progress() as progress:
    #     task_id = progress.add_task(
    #         "Finding lengths from raw FASTA", total=len(sequences)
    #     )
    #     for seq in skbio.read(str(raw_fasta), "fasta", constructor=skbio.DNA):
    #         if seq.metadata["id"] in seq_id_to_new_id:
    #             progress.update(task_id=task_id, advance=1)
    #             sequences[seq_id_to_new_id[seq.metadata["id"]]]["contig_length"] = len(
    #                 seq
    #             )

    #             # include length ratio
    #             sequences[seq_id_to_new_id[seq.metadata["id"]]][
    #                 "contig_length_ratio"
    #             ] = (
    #                 len(seq)
    #                 / sequences[seq_id_to_new_id[seq.metadata["id"]]]["unit_length"]
    #             )

    pd.DataFrame.from_dict(sequences, orient="index")[
        [
            "vdsearch_id",
            "seq_id",
            "description",
            "source",
            "unit_length",
            "gc_content",
            "has_ribozymes",
            "symmetric",
            "rz_plus",
            "rz_minus",
            "rz_plus_from",
            "rz_plus_to",
            "rz_plus_evalue",
            "rz_minus_from",
            "rz_minus_to",
            "rz_minus_evalue",
            "mfe",
            "paired_percent",
            "hairpins",
            "seq",
            "structure",
        ]
    ].to_csv(outfile, sep="\t", index=False, header=header)


@app.callback()
def callback():
    """
    Unstable internal commands.

    ## Beware!

    You are entering the wilderness.
    Beyond the cozy confines of the main pipeline lies mystery, adventure, and, perhaps,
     danger.

    These scripts are not meant to be run by humans, only Benjamin (and maybe Neri).
    It's likely that they are buggy, undocumented, or unrelated to the main pipeline.

    If you're not Benjamin, you should not be here.
    If you're not Neri, you should not be here.
    If you're not either, you should not be here.
    """
    pass


if __name__ == "__main__":
    app()
