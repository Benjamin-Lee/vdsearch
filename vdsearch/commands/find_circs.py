import logging
from pathlib import Path

import nimporter
import typer

from ..nim import find_circs as fc
from ..types import FASTA


def find_circs(
    fasta: Path = FASTA,
    output: Path = typer.Argument(
        ..., file_okay=True, dir_okay=False, writable=True, help="Path to output file"
    ),
):
    """
    Search for circular sequences.

    This method is based on the following paper:

    > Y. Qin *et al.*, “Reference-free and *de novo* Identification of Circular RNAs.”
    > Cold Spring Harbor Laboratory, Apr. 23, 2020. doi: 10.1101/2020.04.21.050617.
    """
    logging.info(f"Searching for circular sequences in {fasta}...")
    count = fc.find_circs(str(fasta), str(output))
    logging.done(f"{count} circRNAs found.")  # type: ignore
