from pathlib import Path

import pandas as pd
import typer

# from rich import print
import rich
from rich.progress import track
from rich.table import Table


def ribozyme_filter(
    infernal_tsv: Path = typer.Argument(
        ...,
        file_okay=True,
        dir_okay=False,
        exists=True,
        readable=True,
        help="Path to Infernal tabular output",
    ),
):

    ribozymes = pd.read_csv(
        infernal_tsv,
        delim_whitespace=True,
        comment="#",
        usecols=[0, 1, 2, 7, 8, 9, 15, 16],
        header=None,
        names=[
            "seq_id",
            "accession",
            "ribozyme",
            "from",
            "to",
            "strand",
            "evalue",
            "inc",
        ],
    )

    double_rz_ids = set(ribozymes.query("strand == '+' & evalue < 0.1").seq_id) & set(
        ribozymes.query("strand == '-' & evalue < 0.1").seq_id
    )
    single_rz_ids = set(ribozymes.query("evalue < 0.01").seq_id) - double_rz_ids
    ribozy_likes_ids = double_rz_ids | single_rz_ids  # union of the two

    # add categorical information about how many ribozymes are in the sequence
    ribozymes.loc[ribozymes["seq_id"].isin(double_rz_ids), "Polarity"] = "(+) and (-)"
    ribozymes.loc[ribozymes["seq_id"].isin(single_rz_ids), "Polarity"] = "(+)"

    single_rzs = ribozymes.loc[ribozymes.seq_id.isin(single_rz_ids)]
    double_rzs = ribozymes.loc[ribozymes.seq_id.isin(double_rz_ids)]
    ribozy_likes = ribozymes.loc[ribozymes.seq_id.isin(ribozy_likes_ids)]

    table = Table(
        title_style="bold magenta",
        title="Ribozyme-containing sequences",
        highlight=True,
    )
    table.add_column("Polarity")
    table.add_column("Count")
    table.add_row("(+)", str(len(single_rz_ids)))
    table.add_row("(+) and (-)", str(len(double_rz_ids)))
    table.add_row("Total", str(len(ribozy_likes_ids)))
    rich.print(table)

    return {
        "single_rzs": single_rzs,
        "double_rzs": double_rzs,
        "ribozy_likes": ribozy_likes,
    }
