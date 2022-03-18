import logging
import os
import subprocess
from pathlib import Path
from pkg_resources import resource_filename

import click
import typer
from vdsearch.types import FASTA
from vdsearch.utils import check_executable_exists, typer_unpacker


@typer_unpacker
def rnamotif(
    fasta: Path = FASTA,
    descr: Path = typer.Argument(..., file_okay=True, dir_okay=False, exists=True),
    output: Path = typer.Argument(
        ...,
        help="Path to TSV output",
        file_okay=True,
        dir_okay=True,
    ),
):
    """
    Run RNAmotif search and format for provided motif descriptions.

    ## Note

    You must have an environment variable `EFNDATA` set to the path of the `efndata` directory from RNAmotif.
    As a convenience (and since the efndata directory hasn't changed in 21 years), we've included it automatically.
    If you provide the path to the efndata directory, it will be used instead.

    ## References

    This method is based on the following paper:

    > Macke, T.J., 2001.
    > RNAMotif, an RNA secondary structure definition and search algorithm.
    > Nucleic Acids Research.
    > https://doi.org/10.1093/nar/29.22.4724
    """

    # TODO: Trigger download of reference CM if not found.

    for tool in ["rnamotif", "rmprune", "rmfmt"]:
        check_executable_exists(
            tool,
        )

    efndata = Path(
        os.environ.get("EFNDATA")
        or resource_filename("vdsearch", "data/rnamotif/efndata")
    )
    logging.debug(f"Using efndata directory: {efndata}")
    if not efndata.exists() or not efndata.is_dir():
        raise FileNotFoundError(
            "Could not find RNAmotif efndata directory. "
            "Please set the environment variable `EFNDATA` to the path of the `efndata` directory from RNAmotif."
        )

    command = (
        f"EFNDATA={efndata} rnamotif "
        f"-descr {descr} {fasta} | "
        f"rmprune | rmfmt > {output}"
    )
    logging.info(f"Searching for {descr.stem} in {fasta} using RNAmotif")
    logging.debug(f"{command=}")
    try:
        subprocess.run(
            command,
            shell=True,
            check=True,
        )
    except subprocess.CalledProcessError as error:
        raise click.Abort(
            f"RNAmotif failed with exit code {error.returncode}",
        )
    logging.done("Done searching for ribozymes.")  # type: ignore
