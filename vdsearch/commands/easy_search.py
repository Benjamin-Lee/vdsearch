import subprocess
import logging
import shutil
from pathlib import Path

import rich_click as click
import pandas as pd
from pkg_resources import resource_filename  # type: ignore
import rich
import skbio
import typer
from rich.console import Console

from vdsearch.commands.canonicalize import canonicalize
from vdsearch.commands.fold import fold
from vdsearch.commands.mmseqs import cluster, search
from vdsearch.commands.dedup import dedup
from vdsearch.commands.find_circs import find_circs
from vdsearch.commands.infernal import infernal
from vdsearch.commands.ribozyme_filter import ribozyme_filter
from vdsearch.commands.rnamotif import rnamotif
from vdsearch.commands.summarize import summarize
from vdsearch.nim import write_seqs as ws
from vdsearch.types import FASTA, ReferenceCms, Threads, ViroidDB
from vdsearch.utils import check_executable_exists


def easy_search(
    fasta: Path = FASTA,
    outdir: Path = typer.Option(
        None,
        file_okay=False,
        dir_okay=True,
        help="Directory to write the results to. If not provided, a directory with the same name as the input file will be created.",
    ),
    reference_db: Path = ViroidDB,
    reference_cms: Path = ReferenceCms,
    circular: bool = typer.Option(
        False, help="Assume circular sequences and skip circularity detection."
    ),
    skip_canonicalization: bool = typer.Option(
        False, help="Skip canonicalization step (useful when pre-annotated)"
    ),
    tmpdir: Path = typer.Option(
        Path("."),
        file_okay=False,
        dir_okay=True,
        help="Path to temporary directory for intermediate files",
    ),
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

    # region: preflight checks
    logging.debug("Checking that all needed tools exist...")
    check_executable_exists("seqkit")
    check_executable_exists("cmsearch")
    check_executable_exists("mmseqs")
    if not reference_db.exists():
        raise click.ClickException(
            f"ViroidDB not found at {reference_db}. Please download it from https://viroids.org."
        )
    if not reference_cms.exists():
        raise click.ClickException(
            f"ReferenceCms not found at {reference_cms}. Please download it using:\n\n\tvdsearch download-cms"
        )

    # handle output directory inference
    if outdir is None:
        logging.debug("No output directory provided. Using input file name.")
        outdir = Path(".") / fasta.stem

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
    # endregion

    # if not reference_db:
    #     download.download_viroiddb()

    # if not reference_cms:
    #     download.download_cms()

    # region: run cirit/rotcanon
    circs = outdir / "circs.fasta"

    if circular and not skip_canonicalization:
        logging.info("Assuming circular sequences.")
        canonicalize(fasta, circs, min_len=100, max_len=10_000)
    elif skip_canonicalization:
        logging.info("Skipping canonicalization step.")
        shutil.copyfile(fasta, circs)

    if not circs.exists():
        find_circs(
            fasta,
            circs,
            canonicalize=True,
            tsv=True,
            min_len=100,
            max_monomer_len=10_000,
            max_len=10_000,
        )
    else:
        logging.warning("CircRNAs already found. Skipping.")
    # endregion

    # region: run dedup using seqkit
    deduped_circs = outdir / "deduped_circs.fasta"
    if not deduped_circs.exists():
        dedup(circs, deduped_circs, threads=threads)
    else:
        logging.warning("CircRNAs already deduplicated. Skipping.")
    # endregion

    # region: run infernal
    cmsearch_output = outdir / "infernal.out"
    cmsearch_tblout = outdir / "infernal.tblout"
    if not cmsearch_output.exists() or not cmsearch_tblout.exists():
        infernal(
            deduped_circs,
            output=cmsearch_output,
            output_tsv=cmsearch_tblout,
            reference_cms=reference_cms,
            threads=threads,
            cmscan=False,
        )
    else:
        logging.warning("Infernal already run. Skipping.")
    # endregion

    # region: run rnamotif if possible
    rnamotif_output = outdir / "rnamotif.tsv"
    if not rnamotif_output.exists():
        try:
            rnamotif(
                deduped_circs,
                Path(resource_filename("vdsearch", "data/rnamotif/Hammerhead_3.descr")),
                rnamotif_output,
            )
        except Exception:
            logging.warn("Could not run RNAmotif. Skipping.")
    else:
        logging.warning("RNAmotif already run. Skipping.")
    # endregion

    # region: find the viroids in the infernal output
    rz_seqs = outdir / "seqs_with_rzs.fasta"
    viroidlike_rzs = outdir / "seqs_with_rzs.tsv"
    if not rz_seqs.exists() or not viroidlike_rzs.exists():
        ribozymes = ribozyme_filter(
            cmsearch_tblout,
            output_tsv=viroidlike_rzs,
            cm_file=reference_cms,
            rnamotif_name="Hammerhead_3",
            rnamotif_txt=rnamotif_output,
        )
        logging.info("Outputting sequences with ribozymes...")
        ws.write_seqs(
            str(deduped_circs),
            str(rz_seqs),
            ribozymes["ribozy_likes"].seq_id.tolist(),
        )
        logging.done(f"Wrote to {rz_seqs}")  # type: ignore
    else:
        logging.warning("Viroid-like sequences already found. Skipping.")
    # endregion

    # region: search against ViroidDB
    viroiddb_hits = outdir / "search_vs_viroiddb.tsv"
    if not viroiddb_hits.exists():
        search(deduped_circs, reference_db, output_tsv=viroiddb_hits, threads=threads)
    # endregion

    # region: extract the viroid matches from ViroidDB
    seqs_matching_viroiddb = outdir / "seqs_matching_viroiddb.fasta"
    if not seqs_matching_viroiddb.exists():
        logging.info("Outputting sequences matching ViroidDB...")
        search_results = pd.read_csv(
            viroiddb_hits,
            sep="\t",
            names="query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,theader,qcov,tcov,cigar".split(
                ","
            ),
        )
        ws.write_seqs(
            str(deduped_circs),
            str(seqs_matching_viroiddb),
            search_results["query"].tolist(),
        )
        logging.done("Wrote to {seqs_matching_viroiddb}")  # type: ignore
    else:
        logging.warning("ViroidDB sequences already emitted. Skipping.")
    # endregion

    # region: merge the two fasta files from ribozyme and viroid db
    merged_seqs = outdir / "viroid_like.fasta"
    if not merged_seqs.exists():
        logging.info("Merging sequences found by ribozyme and ViroidDB searches...")
        seen = set()
        with merged_seqs.open("w") as f:
            # we do viroiddb first so matches are first in the summary
            for search_result in [seqs_matching_viroiddb, rz_seqs]:
                for record in skbio.read(
                    str(search_result), "fasta", constructor=skbio.DNA
                ):
                    if record.metadata["id"] in seen:
                        continue
                    f.write(f">{record.metadata['id']}\n{record}\n")
                    seen.add(record.metadata["id"])
        logging.done("Merged.")  # type: ignore
    else:
        logging.warning("Sequences are already merged. Skipping.")
    # endregion

    # region: fold the sequences
    folded_seqs_plus = outdir / "viroid_like_plus.dbn"
    folded_seqs_minus = outdir / "viroid_like_minus.dbn"
    if not folded_seqs_plus.exists() or not folded_seqs_minus.exists():
        fold(merged_seqs, folded_seqs_plus, threads=threads)

        # reverse complement the sequences and fold them
        rev_comp_temp = outdir / "rev_comp_temp.fasta"
        try:
            subprocess.run(
                f"seqkit seq --reverse --complement --quiet {merged_seqs} --out-file {rev_comp_temp} 2>/dev/null",
                shell=True,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            raise e
        fold(rev_comp_temp, folded_seqs_minus, threads=threads)
        rev_comp_temp.unlink()  # clean up
    else:
        logging.warning("Sequences are already folded. Skipping.")
    # endregion

    # region: generate summary table
    summary_table = outdir / "viroid_like.tsv"
    if not summary_table.exists():
        logging.info("Generating summary table...")
        summary_table = summarize(
            fasta=merged_seqs,
            ribozyme_tsv=viroidlike_rzs,
            viroiddb_tsv=viroiddb_hits,
            dbn_plus=folded_seqs_plus,
            dbn_minus=folded_seqs_minus,
            source=outdir.stem,
            circ_tsv=outdir / "circs.tsv",
            outfile=summary_table,
            header=True,
        )
    else:
        logging.warning("Summary table already generated. Skipping.")
    # endregion

    # region: run clustering
    # mmseqs_tmpdir = tmpdir / "mmseqstmp"
    # cluster_tsv = outdir / "cluster.tsv"
    # rep_seqs = outdir / "cluster.fasta"
    # if not cluster_tsv.exists() or not rep_seqs.exists():
    #     cluster(rz_seqs, tmpdir=mmseqs_tmpdir, threads=threads)
    #     # rename the cluster file to the expected name
    #     # we have to use shutil.move since there might be filesystem issues on clusters
    #     # Python < 3.9 doesn't support Pathlib for shutil.move so we have to use str() :(
    #     shutil.move(str(Path("clu_cluster.tsv")), str(cluster_tsv))
    #     shutil.move(str(Path("clu_rep_seq.fasta")), str(rep_seqs))

    #     # remove the weird pseudo-FASTA formatted cluster file
    #     Path("clu_all_seqs.fasta").unlink()

    #     logging.debug("Cleaning up temporary files from clustering...")
    #     shutil.rmtree(mmseqs_tmpdir)
    # else:
    #     logging.warning("Clustering already run. Skipping.")

    # # report on the clustering results
    # unique_cluster_count = pd.read_csv(
    #     cluster_tsv, sep="\t", names=["cluster_id", "seq_id"]
    # ).cluster_id.nunique()
    # logging.done(f"{unique_cluster_count} viroid-like sequences clusters found.")  # type: ignore
    # endregion

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
