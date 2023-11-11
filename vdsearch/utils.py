import functools
import logging
import shutil
from pathlib import Path
from typing import Callable, Optional

import rich_click as click
import skbio
from typer.models import ParameterInfo


def check_executable_exists(
    name: str, display_name: Optional[str] = None, test_command: str = None
):
    """
    Check if an executable exists in the system path.
    """
    if shutil.which(name) is None:
        logging.error(f"Unable to find {display_name or name} executable.")
        click.rich_click.ERRORS_SUGGESTION = ""
        raise click.ClickException(
            f"Unable to find {display_name or name} executable. Are you sure that it's installed? "
            f"To test if it is, run the following command: {test_command or name}"
        )


def typer_unpacker(f: Callable):
    """
    Make a Typer function into a normal function.

     from https://github.com/tiangolo/typer/issues/279#issuecomment-841875218
    """

    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        default_values = f.__defaults__
        patched_defaults = [
            value.default if isinstance(value, ParameterInfo) else value
            for value in default_values
        ]
        f.__defaults__ = tuple(patched_defaults)

        return f(*args, **kwargs)

    return wrapper

def write_seqs(infile: Path, outfile: Path, seq_ids: list[str]):
    # read each sequence in the input fasta file, write to output fasta file if the sequence id is in the list of seq_ids

    # first convert the list of seq_ids to a set for faster lookup
    seq_ids_set = set(seq_ids)

    # open the input and output files
    with open(infile, "r") as inf, open(outfile, "w") as outf:
        # iterate over the input fasta file using skbio
        for seq in skbio.read(inf, format="fasta", constructor=skbio.DNA):
            # if the sequence id is in the set of seq_ids, write the sequence to the output fasta file
            if seq.metadata["id"] in seq_ids_set:
                skbio.write(seq, format="fasta", into=outf)