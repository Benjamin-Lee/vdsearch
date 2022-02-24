from pathlib import Path
import subprocess
import click
import typer

from vdsearch.types import FASTA, Threads
from ..utils import check_executable_exists
import logging


def dedup(
    fasta: Path = FASTA,
    output: Path = typer.Argument(
        ..., file_okay=True, dir_okay=False, writable=True, help="Path to output file"
    ),
    threads: int = Threads,
):
    """
    Removes duplicate canonicalized sequences.

    Under the hood, this uses seqkit≥2.0.0 and is optimized for canonicalized sequences.
    As such, it is **not safe** to use with non-canonicalized sequences.
    """

    logging.info("Removing duplicate sequences...")

    check_executable_exists(
        "seqkit", display_name="Seqkit", test_command="seqkit version"
    )

    # TODO: check seqkit version is >= 2.0.0

    command = (
        "seqkit rmdup "
        "--by-seq "
        "--only-positive-strand "
        f"--threads {threads} "
        f"--ignore-case {fasta} -o {output}"
    )
    logging.debug("Running command [cyan]{}[/]".format(command))

    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True)
        count = result.stderr.decode("utf-8").split()[1]
        logging.done(f"{count} duplicate sequences removed.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed with exit code {e.returncode}")
        logging.error(f"Command output: '{e.output.decode('utf-8').rstrip()}'")
        raise click.Abort("Command failed. See logs for error messages.")
