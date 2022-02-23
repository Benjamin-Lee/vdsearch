import typer

FASTA = typer.Argument(
    ...,
    help="Path to .fasta or .fasta.gz file",
    exists=True,
    file_okay=True,
    dir_okay=False,
    readable=True,
)
