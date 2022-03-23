import os
from pathlib import Path
import typer

FASTA = typer.Argument(
    ...,
    help="Path to .fasta file",
    exists=True,
    file_okay=True,
    dir_okay=False,
    readable=True,
)
Threads: int = typer.Option(
    os.cpu_count(),
    help="Threads to use when applicable. Default is all available.",
)

ReferenceCms: Path = typer.Option(
    Path(typer.get_app_dir("vdsearch")) / "data" / "cms" / "merged",
    help="Path to Infernal ribozyme families. If none is provided, the latest ViroidDB will be used.",
    file_okay=True,
    dir_okay=False,
)
