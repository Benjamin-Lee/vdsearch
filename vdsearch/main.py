from typing import List

import typer
from halo import Halo
import nimporter
import subprocess
from yaspin import yaspin

from .internal import internal

from .rich_wrapper import MyTyper


# nimporter.IGNORE_CACHE = True  # Rebuild all Nim files

from vdsearch.find_circs import find_circs as fc


import time

from pathlib import Path

app = MyTyper()

FASTA = typer.Argument(
    ...,
    help="Path to .fasta or .fasta.gz file",
    exists=True,
    file_okay=True,
    dir_okay=False,
    readable=True,
)


@app.command()
def download_viroiddb():
    """Download the latest ViroidDB dataset."""

    with Halo(text="Downloading ViroidDB dataset...") as sp:
        time.sleep(2)
        sp.succeed("Downloaded ViroidDB")


@app.command()
def download_cms():
    """Download the latest covariance matrices for ribozymes."""

    with Halo(text="Downloading ribozyme CMs dataset...") as sp:
        time.sleep(2)
        sp.succeed("Downloaded ribozyme CMs")


@app.command()
def find_circs(
    fasta: Path = FASTA,
    output: Path = typer.Argument(
        ..., file_okay=False, dir_okay=False, writable=True, help="Path to output file"
    ),
):
    """Find circular RNA sequences."""
    with Halo(text=f"Finding circular regions in {fasta}...") as sp:
        time.sleep(2)
        # fc(str(fasta.absolute()), str(output.absolute()))
        sp.succeed("Found circular contigs")
    pass


@app.command()
def dedup(fasta: Path = FASTA):
    """Remove duplicate sequences."""
    with Halo(text=f"Removing duplicate sequences from {fasta}...") as sp:
        time.sleep(2)
        sp.succeed("Removed duplicate sequences")
    pass


@app.command()
def infernal(fasta: Path = FASTA):
    """Run Infernal on the sequences."""
    with Halo(text=f"Running Infernal on {fasta}...") as sp:
        time.sleep(2)
        sp.succeed("Ran Infernal")
    pass


@app.command()
def mmseqs(fasta: Path = FASTA):
    """Run MMseqs on the sequences."""
    with Halo(text=f"Running MMseqs on {fasta}...") as sp:
        time.sleep(2)
        sp.succeed("Ran MMseqs")
    pass


@app.command()
def remove_rz_only_hits(fasta: Path = FASTA):
    """Remove sequences that only hit ribozymes."""
    with Halo(text=f"Removing sequences that only hit ribozymes from {fasta}...") as sp:
        time.sleep(2)
        sp.succeed("Generated new FASTA file of viroid-like sequences")
    pass


@app.command()
def cluster(fasta: Path = FASTA):
    """Cluster circular sequences."""
    logfile = Path("mmseqs.log.txt")
    typer.secho("ℹ ", fg="cyan", nl=False)
    typer.echo(f"Logs are in {logfile}")
    with yaspin(text=f"Clustering {fasta.name}", color="cyan", timer=True) as sp:
        subprocess.run(
            # f"seqkit concat --quiet {fasta} {fasta} | "
            " mmseqs easy-cluster "
            "--min-seq-id 0.75  "
            "--min-aln-len 100 "
            "--seq-id-mode 2 "
            "--cov-mode 0 "
            "-e 0.0001 "
            "-c 0.75 "
            "-s 7.5 "
            "--cluster-mode 1 "
            f"{fasta} clu tmp > {logfile}",
            shell=True,
        )
        sp.color = "green"
        sp.text = f"Clustered {fasta.name}"
        sp.ok("✔")
    return 1


@app.command()
def deconcat(fasta: Path = FASTA):
    """Deconcatenate circular sequences.

    Sometimes, it makes sense to concatenate circular sequences to themselves for processing using a tool that assumes linearity.

    """
    with Halo(text=f"Deconcatenating {fasta}...") as sp:
        sp.succeed("Deconcatenated")
    pass


@app.command()
def easy_search(
    fasta: Path = FASTA,
    reference_db: Path = typer.Option(
        None,
        help="Path to FASTA-formatted reference viroid database. If none is provided, the latest ViroidDB will be used.",
        file_okay=True,
        dir_okay=False,
    ),
    reference_cms: Path = typer.Option(
        None,
        help="Path to Infernal ribozyme families. If none is provided, the latest ViroidDB will be used.",
        file_okay=True,
        dir_okay=False,
    ),
):
    """
    Search for viroid-like sequences.
    """

    if not reference_db:
        download_viroiddb()

    if not reference_cms:
        download_cms()

    # run cirit/rotcanon/dedupe
    find_circs(fasta)

    dedup(fasta)

    # run infernal
    infernal(fasta)

    # run mmseqs
    mmseqs(fasta)

    remove_rz_only_hits(fasta)

    typer.secho("\nDone!", fg="green", nl=False)
    typer.secho(
        f" The results are in {fasta.stem}_viroidlikes.fasta.",
    )


# Support for using --version
# See https://github.com/tiangolo/typer/issues/52
__version__ = "0.0.1"


def version_callback(value: bool):
    if value:
        typer.echo(f"vdsearch {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Show the version and exit.",
    )
):
    """
    Workflows and utilities for finding and analyzing viroid-like circular RNAs.
    """
    return


app.add_typer(internal.app, name="internal")
if __name__ == "__main__":
    app(prog_name="vdsearch")
