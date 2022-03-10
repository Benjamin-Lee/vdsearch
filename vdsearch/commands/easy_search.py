import logging
import shutil
from pathlib import Path

import click
import pandas as pd
import rich
import typer
from rich.console import Console

from vdsearch.commands.canonicalize import canonicalize
from vdsearch.commands.cluster import cluster
from vdsearch.commands.dedup import dedup
from vdsearch.commands.find_circs import find_circs
from vdsearch.commands.infernal import infernal
from vdsearch.commands.ribozyme_filter import ribozyme_filter
from vdsearch.nim import write_seqs as ws
from vdsearch.types import FASTA, ReferenceCms, Threads
from vdsearch.utils import check_executable_exists


def easy_search(
    fasta: Path = FASTA,
    outdir: Path = typer.Option(Path("."), file_okay=False, dir_okay=True),
    reference_db: Path = typer.Option(
        None,
        help="Path to FASTA-formatted reference viroid database. If none is provided, the latest ViroidDB will be used.",
        file_okay=True,
        dir_okay=False,
    ),
    reference_cms: Path = ReferenceCms,
    assume_circular: bool = typer.Option(False, help="Assume circular sequences"),
    threads: int = Threads,
    force: bool = typer.Option(
        False,
        help="Force running even if lockfile is present. Internal and inadvisable for production.",
        hidden=True,
    ),
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

    # preflight checks
    logging.debug("Checking that all needed tools exist...")
    check_executable_exists("seqkit")
    check_executable_exists("cmsearch")
    check_executable_exists("mmseqs")

    logging.debug("Making output directory...")
    outdir.mkdir(parents=True, exist_ok=True)

    # add a lock to the outdir to prevent parallel runs or incomplete data.
    lockfile = outdir / ".vdsearch.lock"
    if lockfile.exists() and not force:
        raise click.ClickException(
            f"Lockfile {lockfile} exists in the output directory.\n\n"
            "Probably, a previous run is still in progress or has failed. "
            f"If it has failed, data in {outdir} may be corrupt or incomplete. "
            "Please wait for the previous run to finish or delete the lockfile if you are sure the data is not corrupt.\n\n"
            f"To delete the lockfile, run:\n\n\trm {lockfile}"
        )
    lockfile.touch()
    logging.debug("Directory locked.")

    logging.info(f"Beginning search for viroid-like RNAs using {threads} threads...")
    # if not reference_db:
    #     download.download_viroiddb()

    # if not reference_cms:
    #     download.download_cms()

    # run cirit/rotcanon
    circs = outdir / Path(f"01.{fasta.stem}.circs.fasta")

    if assume_circular:
        logging.info("Assuming circular sequences.")
        canonicalize(fasta, circs)

    if not circs.exists():
        find_circs(fasta, circs, canonicalize=True)
    else:
        logging.warning("CircRNAs already found. Skipping.")

    # run dedup using seqkit
    deduped_circs = outdir / Path(f"03.{fasta.stem}.deduped.fasta")
    if not deduped_circs.exists():
        dedup(circs, deduped_circs, threads=threads)
    else:
        logging.warning("CircRNAs already deduplicated. Skipping.")

    # run infernal
    cmsearch_output = outdir / Path(f"04.{fasta.stem}.infernal.out")
    cmsearch_tsv = outdir / Path(f"04.{fasta.stem}.infernal.tsv")
    if not cmsearch_output.exists() and not cmsearch_tsv.exists():
        infernal(
            deduped_circs,
            output=cmsearch_output,
            output_tsv=cmsearch_tsv,
            reference_cms=reference_cms,
            threads=threads,
            cmscan=False,
        )
    else:
        logging.warning("Infernal already run. Skipping.")

    # find the viroids in the infernal output
    rz_seqs = outdir / Path(f"05.{fasta.stem}.ribozymes.fasta")
    if not rz_seqs.exists():
        ribozymes = ribozyme_filter(cmsearch_tsv, cm_file=reference_cms)
        logging.info("Outputting viroid-like sequences...")
        ws.write_seqs(
            str(deduped_circs),
            str(rz_seqs),
            ribozymes["ribozy_likes"].seq_id.tolist(),
        )
        logging.done(f"Wrote to {rz_seqs}")  # type: ignore
    else:
        logging.warning("Viroid-like sequences already found. Skipping.")

    # run clustering
    mmseqs_tmpdir = Path("mmseqstmp")
    cluster_tsv = outdir / Path(f"06.{fasta.stem}.cluster.tsv")
    if not cluster_tsv.exists():
        cluster(rz_seqs, tmpdir=mmseqs_tmpdir, threads=threads)
        # rename the cluster file to the expected name
        # we have to use shutil.move since there might be filesystem issues on clusters
        # Python < 3.9 doesn't support Pathlib for shutil.move so we have to use str() :(
        shutil.move(str(Path("clu_cluster.tsv")), str(cluster_tsv))
        logging.debug("Cleaning up temporary files from clustering...")
        shutil.rmtree(mmseqs_tmpdir)
    else:
        logging.warning("Clustering already run. Skipping.")

    # report on the clustering results
    unique_cluster_count = pd.read_csv(
        cluster_tsv, sep="\t", names=["cluster_id", "seq_id"]
    ).cluster_id.nunique()
    logging.done(f"{unique_cluster_count} viroid-like sequences clusters found.")  # type: ignore

    # # run mmseqs
    # mmseqs(fasta)

    # remove_rz_only_hits(fasta)

    # logging.done(f"The results are in [green bold]{fasta.stem}_viroidlikes.fasta.[/]")  # type: ignore
    Console().log(
        # "\n",
        # "Thanks for using [green]vdsearch[/]!",
        "\n",
        rich.panel.Panel(
            rich.markdown.Markdown(
                """
If you use these results in your research, please cite:

> B.D. Lee *et al.* (2022) vdsearch: A tool for viroid-like RNA searches.
            """
            ),
            title_align="left",
            border_style="dim",
            width=88,
            title="Citation",
        ),
        "\n",
        "[dim]Brought to you by: [cyan bold]NIH/NLM/NCBI[/], [blue bold]University of Oxford[/], and [bold]Tel Aviv University[/][/dim]",
        "\n",
    )

    # remove the lock
    lockfile.unlink()
    logging.debug("Lockfile removed.")
