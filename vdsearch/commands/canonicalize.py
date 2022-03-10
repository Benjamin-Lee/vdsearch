from pathlib import Path
from rich.console import Console

import typer

from vdsearch.types import FASTA
from vdsearch.nim import canonicalize as rotcanon
from vdsearch.utils import typer_unpacker


@typer_unpacker
def canonicalize(
    fasta: Path = FASTA,
    output: Path = typer.Argument(
        ..., file_okay=True, dir_okay=False, writable=True, help="Path to output file"
    ),
):
    """Compute a rotationally canonical representation of each sequence in a FASTA file."""
    console = Console()
    with console.status("Computing canonical representations..."):
        rotcanon.canonicalize(str(fasta), str(output))
