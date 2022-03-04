import logging
from pathlib import Path
from typing import Optional

import nimporter
import pandas as pd
import typer
from rich.console import Console
from vdsearch.nim import write_seqs as ws
from vdsearch.rich_wrapper import MyTyper
from vdsearch.types import FASTA

app = MyTyper(hidden=True)


@app.command()
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
                        [[seq_id, structure, seq, mfe, paired_percent, hairpins]],
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
