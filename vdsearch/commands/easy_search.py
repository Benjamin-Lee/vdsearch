import logging
from pathlib import Path
from turtle import title
import rich
from rich.console import Console

import typer
from vdsearch.commands.canonicalize import canonicalize

from vdsearch.types import FASTA, ReferenceCms, Threads
from . import (
    download_cms,
    download_viroiddb,
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
    reference_cms: Path = ReferenceCms,
    threads: int = Threads,
):
    """
    Search for viroid-like sequences.


    ## Pipeline description

    At a high level, this workflow follows these steps:

    1. Find circRNAs
    2. For each circRNA compute a canonical representation
    3. Deduplicate the circRNAs
    4. Search them for ribozymes
    5. Also search the circRNAs against a database of known viroid-like RNAs (ViroidDB)
    6. Using the ribozyme data and search results, output viroid-like sequences
    """

    logging.info("Beginning search for viroid-like RNAs...")
    # if not reference_db:
    #     download.download_viroiddb()

    # if not reference_cms:
    #     download.download_cms()

    # run cirit/rotcanon/dedupe
    circs = f"01.{fasta.stem}.circs.fasta"
    find_circs(Path(fasta), Path(circs))

    canon = f"02.{fasta.stem}.canon.fasta"
    canonicalize(Path(circs), Path(canon))

    deduped_circs = f"03.{fasta.stem}.deduped.fasta"
    dedup(Path(canon), Path(deduped_circs), threads=threads)

    # run infernal
    cmsearch_output = f"04.{fasta.stem}.infernal.out"
    cmsearch_tsv = f"04.{fasta.stem}.infernal.tsv"
    infernal(
        Path(deduped_circs),
        Path(cmsearch_output),
        Path(cmsearch_tsv),
        reference_cms=reference_cms,
        threads=threads,
    )

    # # run mmseqs
    # mmseqs(fasta)

    # remove_rz_only_hits(fasta)

    logging.done(f"The results are in [green bold]{fasta.stem}_viroidlikes.fasta.[/]")
    Console().print(
        "\n",
        rich.panel.Panel(
            rich.markdown.Markdown(
                """
If you use these results in your research, please cite:

> B.D. Lee et al. (2022) vdsearch: A tool for viroid-like RNA searches.
            """
            ),
            title_align="left",
            border_style="dim",
            width=88,
            title="Citation",
        ),
    )
    return "done"
