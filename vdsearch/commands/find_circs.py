import logging
from pathlib import Path
import sys

import nimporter
import typer

from vdsearch.nim import find_circs as fc


def find_circs(
    fasta: Path = typer.Argument(
        ...,
        help="Path to .fasta file or '-' for stdin",
        file_okay=True,
        dir_okay=False,
        readable=True,
    ),
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
        1,
        help="Minimum length of the sequence to be monomerized and included in the output",
    ),
    max_len: int = typer.Option(
        sys.maxsize,
        help="Maximum length of the input sequence to be monomerized and included in the output",
    ),
    max_monomer_len: int = typer.Option(
        sys.maxsize,
        help="Maximum length of the resulting monomer to be included in the output",
    ),
):
    """
    Search for circular sequences.

    ## Performance notes

    By default, the full search routine is run for every sequence.
    This can be slow for large datasets where many of the huge input sequences couldn't possibly be circRNAs.
    To speed up the search, you can provide a limit on the smallest (**--min-len**) and largest input sequence lengths (**--max-len**) to consider and report.
    For reporting, you can also set the maximum reported monomer length with **--max-monomer-len**.

    ## References

    This method is based on the following paper:

    > Qin, Y., Xu, T., Lin, W., Jia, Q., He, Q., Liu, K., et al., 2020.
    > Reference-free and *de novo* Identification of Circular RNAs (preprint).
    > <https://doi.org/10.1101/2020.04.21.050617>

    """
    logging.info(f"Searching for circular sequences in {fasta}.")
    logging.debug(f"{'Will' if canonicalize else 'Will not'} canonicalize them.")
    logging.debug(
        f"Only monomers of length {min_len} to {max_len} will be monomerized. "
        f"Only monomers less than {max_monomer_len} nt will be d."
    )
    count = fc.find_circs(
        str(fasta),
        str(output),
        canonicalize=canonicalize,
        outTsv=tsv,
        minLen=min_len,
        maxLen=max_len,
        maxMonomerLen=max_monomer_len,
        verbose=logging.getLogger().isEnabledFor(logging.DEBUG),
    )
    logging.done(f"{count[0]:,} circRNAs found from {count[1]:,} sequences ({(count[2] / 1000000) / (count[3] / 1000.0):,.0f} Mbp/sec).")  # type: ignore
