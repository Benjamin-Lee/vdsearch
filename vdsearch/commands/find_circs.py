import logging
from pathlib import Path
import sys

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
    tsv: bool = typer.Option(
        True,
        help="Output a .tsv file with the ratio between the length of the original sequence and its monomerized version. If present, the output file will be named <output>.tsv",
    ),
    min_len: int = typer.Option(
        1, help="Minimum length of the monomer to be included in the output"
    ),
    max_len: int = typer.Option(
        sys.maxsize, help="Maximum length of the monomer to be included in the output"
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
    logging.debug(f"Only monomers of length {min_len} to {max_len} will be included.")
    count = fc.find_circs(
        str(fasta),
        str(output),
        canonicalize=canonicalize,
        outTsv=tsv,
        minLen=min_len,
        maxLen=max_len,
    )
    logging.done(f"{count[0]:,} circRNAs found from {count[1]:,} sequences ({count[1] / (count[2] / 1000.0):,.0f} seqs/sec).")  # type: ignore
