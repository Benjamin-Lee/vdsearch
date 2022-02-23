from pathlib import Path

import typer

from vdsearch.types import FASTA
from . import (
    download,
    find_circs,
    dedup,
    cluster,
    infernal,
)


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
        download.download_viroiddb()

    if not reference_cms:
        download.download_cms()

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
