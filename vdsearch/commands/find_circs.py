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

    > Qin, Yangmei, Tingting Xu, Wenbo Lin, Qingjie Jia, Qiushun He, Ke Liu, Juan Du, *et al*.
    > “Reference-Free and *de Novo* Identification of Circular RNAs”
    > (Cold Spring Harbor Laboratory, April 2020).
    > [https://doi.org/10.1101/2020.04.21.050617](https://doi.org/10.1101/2020.04.21.050617).

    """
    logging.info(f"Searching for circular sequences in {fasta}.")
    logging.debug(f" {'Will' if canonicalize else 'Will not'} canonicalize them.")
    count = fc.find_circs(str(fasta), str(output), canonicalize=canonicalize)
    logging.done(f"{count} circRNAs found.")  # type: ignore
