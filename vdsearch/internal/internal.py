import hashlib
import logging
from pathlib import Path
from typing import List, Optional

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
def dbn2dot(
    dbn: Path,
    outfile: Path = typer.Argument(..., help="Path to output file"),
):
    """
    Convert a dbn file to a dot file
    """
    seq_id, structure, seq, mfe = "", "", "", 0.0
    with dbn.open() as f:
        for i, line in enumerate(f):
            if i % 3 == 0:
                seq_id = line.strip()[1:]
            elif i % 3 == 1:
                seq = line.strip()
            elif i % 3 == 2:
                structure = line.strip()[0 : line.strip().rfind("(")]
            if i == 3:
                break

    # print("digraph G {\n\t{")
    # for i, base in enumerate(seq):
    #     print(f"\t\t{i} [label={base}]")
    # print("\t}")
    # for i in range(len(seq) - 1):
    #     print(f"{i} -> {i + 1}")
    # print(f"{len(seq) - 1} -> 0")

    open = []
    for i, base in enumerate(structure):
        if base == "(":
            open.append(i)
        elif base == ")":
            # print(f"{open.pop()} -> {i}[dir=none]")
            print(f"{open.pop()}\t{i}")

    # print("}")


@app.command()
def merge_summaries(
    infiles: List[Path] = typer.Argument(
        ...,
        help="Path to summary files. If a directory, the default summary file from easy-search is used.",
    ),
    outfile: Path = typer.Argument(..., help="Path to output file"),
):
    """
    Merge multiple summary files into a single file.
    """
    # read in the files
    dfs = []
    for infile in infiles:
        # allow the user to specify a directory or a file
        # if it's a dir, use the default file name from easy-search
        if infile.is_dir():
            infile = infile / "viroid_like.tsv"
        dfs.append(pd.read_csv(infile, sep="\t"))

    # merge the dataframes
    df = pd.concat(dfs)

    # write the output file
    df.to_csv(outfile, sep="\t", index=False)


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
