import logging
import subprocess
from pathlib import Path

import click
import typer

from vdsearch.types import FASTA, Threads
from vdsearch.utils import check_executable_exists, typer_unpacker


@typer_unpacker
def cluster(
    fasta: Path = FASTA,
    tmpdir: Path = typer.Option(
        Path("tmp"),
        help="Path to temporary directory to use for intermediate files",
    ),
    threads=Threads,
):
    """Cluster circular sequences."""
    check_executable_exists("mmseqs")

    logfile = Path("mmseqs.log.txt")

    logging.info("Clustering...")
    command = (
        # f"seqkit concat --quiet {fasta} {fasta} | "
        "mmseqs easy-cluster "
        "--min-seq-id 0.75 "
        "--min-aln-len 100 "
        "--seq-id-mode 2 "
        "--cov-mode 0 "
        "-e 0.0001 "
        "-c 0.75 "
        "-s 7.5 "
        "--cluster-mode 1 "
        f"--threads {threads} "
        f"{fasta} clu {tmpdir} "
        f"> {logfile}"
    )
    logging.debug(f"{command=}")
    try:
        subprocess.run(
            command,
            shell=True,
            check=True,
        )
        logging.done("Done clustering.")  # type: ignore
    except subprocess.CalledProcessError as e:
        logging.error(f"MMseqs failed with exit code {e.returncode}")
        if e.output:
            logging.error(f"MMseqs output: '{e.output.decode('utf-8').rstrip()}'")
        raise click.ClickException("MMseqs failed. See logs for error messages.")


@typer_unpacker
def search(
    query: Path = FASTA,
    target: Path = FASTA,
    output_tsv: Path = typer.Argument(..., file_okay=True, dir_okay=False),
    tmpdir: Path = typer.Option(
        Path("tmp"),
        help="Path to temporary directory to use for intermediate files",
    ),
    threads=Threads,
):
    """Search sequences using MMseqs.

    There's nothing fancy going on here, just a wrapper around MMseqs.
    """
    check_executable_exists("mmseqs")

    logfile = Path("mmseqs.log.txt")

    logging.info(f"Searching against {target.name}...")
    command = (
        "mmseqs easy-search "
        "-s 7.5 "
        "--search-type 3 "
        "--format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,theader,qcov,tcov,cigar "
        f"--threads {threads} "
        f"{query} '{target}' {output_tsv} {tmpdir} "
        f"> {logfile}"
    )
    logging.debug(f"{command=}")
    try:
        subprocess.run(
            command,
            shell=True,
            check=True,
        )
        logging.done("Done searching..")  # type: ignore
    except subprocess.CalledProcessError as e:
        logging.error(f"MMseqs failed with exit code {e.returncode}")
        if e.output:
            logging.error(f"MMseqs output: '{e.output.decode('utf-8').rstrip()}'")
        raise click.ClickException("MMseqs failed. See logs for error messages.")
