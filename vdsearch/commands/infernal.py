import logging
import subprocess
from pathlib import Path

import click

from ..utils import check_executable_exists
import typer
from rich.console import Console
from vdsearch.types import FASTA, ReferenceCms, Threads


def infernal(
    fasta: Path = FASTA,
    output: Path = typer.Argument(
        ...,
        help="Path to Infernal text output",
        file_okay=True,
        dir_okay=False,
    ),
    output_tsv: Path = typer.Argument(
        ...,
        help="Path to parseable Infernal tabular output",
        file_okay=True,
        dir_okay=False,
    ),
    reference_cms: Path = ReferenceCms,
    threads: int = Threads,
):
    """
    Run Infernal cmsearch for provided covariance matricies.
    """
    if reference_cms is None:
        raise click.Abort(
            "No reference CM provided. "
            "Currently unable to automatically retrieve one."
        )
    check_executable_exists("cmsearch", test_command="cmsearch -h")
    console = Console()
    with console.status("Searching using Infernal..."):
        subprocess.run(
            f"cmsearch --cpu {threads} --tblout {output_tsv} -o {output} {reference_cms} {fasta}",
            shell=True,
            check=True,
        )
