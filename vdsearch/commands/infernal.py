import logging
import subprocess
from pathlib import Path
from typing import Optional

import click
import typer
from vdsearch.types import FASTA, ReferenceCms, Threads
from vdsearch.utils import check_executable_exists, typer_unpacker


@typer_unpacker
def infernal(
    fasta: Path = FASTA,
    reference_cms: Path = ReferenceCms,
    output: Optional[Path] = typer.Option(
        None,
        help="Path to Infernal text output",
        file_okay=True,
        dir_okay=False,
    ),
    output_tsv: Optional[Path] = typer.Option(
        None,
        help="Path to parseable Infernal tabular output.",
        file_okay=True,
        dir_okay=False,
    ),
    cmscan: bool = typer.Option(
        False, "--cmscan/--cmsearch", help="Infernal program to use"
    ),
    cut_ga: bool = typer.Option(
        False,
        help="Use the GA (gathering) bit scores in the model to set hit reporting and inclusion thresholds. GA thresholds are generally considered to be the reliable curated thresholds defining family membership; for example, in Rfam, these thresholds define what gets included in Rfam Full alignments based on searches with Rfam Seed models.",
    ),
    cut_nc: bool = typer.Option(
        False,
        help="Use the NC (noise cutoff) bit score thresholds in the model to set hit reporting and inclusion thresholds. NC thresholds are generally considered to be the score of the highest-scoring known false positive.",
    ),
    cut_tc: bool = typer.Option(
        False,
        help="Use the TC (trusted cutoff) bit score thresholds in the model to set hit reporting and inclusion thresholds. TC thresholds are generally considered to be the score of the lowest-scoring known true positive that is above all known false positives.",
    ),
    evalue: Optional[float] = typer.Option(
        None, help="The maximum E value to report.", min=0.0
    ),
    threads: int = Threads,
):
    """
    Run Infernal cmsearch or cmscan for provided covariance matrices.

    ## Note

    Descriptions for **--cut-nc**, **--cut-tc** and **--cut-ga** are copied directly from [Infernal's manpage](http://eddylab.org/infernal/Userguide.pdf).
    Note also that, unlike in normal Infernal commands, the flags are kebab-case rather than snake_case.

    ## References

    This method is based on the following paper:

    > Nawrocki, E. P., and S. R. Eddy.
    > “Infernal 1.1: 100-Fold Faster RNA Homology Searches.”
    > Bioinformatics (Oxford University Press (OUP), September 2013).
    > [https://doi.org/10.1093/bioinformatics/btt509](https://doi.org/10.1093/bioinformatics/btt509).
    """

    base_command = "cmscan" if cmscan else "cmsearch"

    if reference_cms is None:
        raise click.Abort(
            "No reference CM provided. "
            "Currently unable to automatically retrieve one."
        )

    # TODO: Trigger download of reference CM if not found.

    check_executable_exists(
        base_command,
        test_command=f"{base_command} -h",
    )

    # TODO: check if input is a valid FASTA file

    # default output files if none provided
    if output is None:
        output = Path(fasta.stem + "." + base_command + ".out")
    if output_tsv is None:
        output_tsv = Path(fasta.stem + "." + base_command + ".tsv")

    command = (
        f"{base_command} "
        f"--cpu {threads} "
        f"--tblout {output_tsv} "
        f"{'--cut_ga ' if cut_ga else ''}"
        f"{'--cut_tc ' if cut_tc else ''}"
        f"{'--cut_nc ' if cut_nc else ''}"
        f"{'-E ' + str(evalue) + ' ' if evalue is not None and float(evalue) else ''}"
        f"-o {output} '{reference_cms}' {fasta}"
    )
    logging.info(f"Searching {fasta} using {base_command}")
    logging.debug(f"{command=}")
    try:
        subprocess.run(
            command,
            shell=True,
            check=True,
        )
    except subprocess.CalledProcessError as error:
        raise click.Abort(
            f"{base_command} failed with exit code {error.returncode}",
        )
    logging.done("Done searching for ribozymes.")  # type: ignore
