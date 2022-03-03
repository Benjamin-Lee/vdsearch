import logging
from pathlib import Path

import nimporter
import typer

from vdsearch.nim import find_circs as fc
from vdsearch.types import FASTA


def find_circs(
    fasta: Path = FASTA,
    output: Path = typer.Argument(
        ..., file_okay=True, dir_okay=False, writable=True, help="Path to output file"
    ),
    canonicalize: bool = typer.Option(
        True, help="Output canonicalized sequences instead of raw sequences"
    ),
):
    """
    Search for circular sequences.

    ## References

    This method is based on the following paper:

    > Qin, Y., Xu, T., Lin, W., Jia, Q., He, Q., Liu, K., et al., 2020.
    > Reference-free and *de novo* Identification of Circular RNAs (preprint).
    > <https://doi.org/10.1101/2020.04.21.050617>

    """
    logging.info(f"Searching for circular sequences in {fasta}.")
    logging.debug(f"{'Will' if canonicalize else 'Will not'} canonicalize them.")
    count = fc.find_circs(str(fasta), str(output), canonicalize=canonicalize)
    logging.done(f"{count} circRNAs found.")  # type: ignore
