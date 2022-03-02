from pathlib import Path
import nimporter
import pandas as pd
import typer
from vdsearch.types import FASTA
from vdsearch.nim import write_seqs as ws

from vdsearch.utils import check_executable_exists


def write_seqs(
    viroidlike_tsv: Path = typer.Argument(
        ...,
        help="Path to parsed Infernal .tsv file from ribozyme-search",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
    ),
    input_fasta: Path = FASTA,
    output: Path = typer.Argument(..., file_okay=True, dir_okay=False, writable=True),
):
    """
    Write viroid-like sequences identified using ribozyme-filter.

    You probably don't need to use this command directly.
    It should only be used if you are calling ribozymes-filter directly.
    """

    ids = (
        pd.read_csv(viroidlike_tsv, sep="\t", usecols=["seq_id"])["seq_id"]
        .unique()
        .tolist()
    )

    ws.write_seqs(str(input_fasta), str(output), ids)
