from pathlib import Path
import subprocess
from sys import stderr
from rich.console import Console
import typer

from vdsearch.types import FASTA
import logging


def dedup(
    fasta: Path = FASTA,
    output: Path = typer.Argument(
        ..., file_okay=True, dir_okay=False, writable=True, help="Path to output file"
    ),
):
    """
    Removes duplicate canonicalized sequences.

    Under the hood, this uses seqkit > 2.0.0.
    """

    console = Console()
    logging.warn(
        "This command expects canonicalized and monomerized sequences.",
    )
    with console.status("Removing duplicates..."):
        command = (
            "seqkit rmdup "
            "--by-seq "
            "--only-positive-strand "
            f"--ignore-case {fasta} -o {output}"
        )
        logging.debug("Running command [cyan]{}[/]".format(command))
        try:
            result = subprocess.run(
                command, shell=True, check=True, capture_output=True
            )
            count = result.stderr.decode("utf-8").split()[1]
            console.log(
                "[green]:heavy_check_mark:[/] Done removing duplicates. "
                f"{count} sequences removed."
            )
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed with exit code {e.returncode}")
            logging.error(f"Command output: '{e.output.decode('utf-8').rstrip()}'")
