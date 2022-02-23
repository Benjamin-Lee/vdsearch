import logging
from pathlib import Path
import subprocess
from numpy import absolute
from rich.console import Console

import typer

from vdsearch.types import FASTA


def cluster(fasta: Path = FASTA):
    """Cluster circular sequences."""
    logfile = Path("mmseqs.log.txt")
    console = Console()
    with console.status(f"Clustering [cyan]{str(fasta)}[/]"):
        try:
            subprocess.run(
                f"seqkit concat --quiet {fasta} {fasta} | "
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
                check=True,
            )
            console.log("[green]:heavy_check_mark:[/] Done clustering.")
        except subprocess.CalledProcessError as e:
            logging.error(f"Command failed with exit code {e.returncode}")
            logging.error(f"Command output: '{e.output.decode('utf-8').rstrip()}'")
