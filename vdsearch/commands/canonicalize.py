import logging
from pathlib import Path
import sys
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
    min_len: int = typer.Option(
        1,
        help="Minimum length of a sequence to be canonicalized and included in the output",
    ),
    max_len: int = typer.Option(
        sys.maxsize,
        help="Maximum length of a sequence to be canonicalized and included in the output",
    ),
):
    """Compute a rotationally canonical representation of each sequence in a FASTA file.

    ## Performance notes

    By default, every sequence is canonicalized.
    This can be slow for large sequences since canonicalization is O(n^2).
    To speed things up, you can provide a limit on the smallest (**--min-len**) and largest input sequence lengths (**--max-len**) to output.
    """
    logging.info(f"Canonicalizing sequences in {fasta}.")
    logging.debug(
        f"Only sequences of length {min_len} to {max_len} will be canonicalized. "
    )
    rotcanon.canonicalize(str(fasta), str(output), minLen=min_len, maxLen=max_len)
    logging.done("Done canonicalizing.")  # type: ignore
