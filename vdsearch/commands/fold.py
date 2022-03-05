import logging
from pathlib import Path
import subprocess
from typing import Optional

import typer

from vdsearch.types import FASTA, Threads
from vdsearch.utils import check_executable_exists


def fold(
    fasta: Path = FASTA,
    output: Optional[Path] = typer.Option(None),
    threads: int = Threads,
    temp: int = typer.Option(25, help="Number of degrees C to use for folding"),
):
    """Predict secondary structures of circRNAs

    ## References

    > Lorenz, R., Bernhart, S.H., Höner zu Siederdissen, C., Tafer, H., Flamm, C., Stadler, P.F., et al., 2011.
    > ViennaRNA Package 2.0.
    > Algorithms Mol Biol 6, 26.
    > https://doi.org/10.1186/1748-7188-6-26

    """
    check_executable_exists("rnafold", "ViennaRNA")

    if output is None:
        logging.debug("Using default output filename.")
        output = Path(fasta).with_suffix(".dbn")

    command = f"rnafold --circ --noPS --jobs={threads} --temp={temp} {fasta} > {output}"
    logging.debug(f"{command=}")
    logging.info(f"Folding {fasta}")
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(
            f"Command [cyan]{command}[/] failed with exit code [red]{e.returncode}[/]"
        )
        raise e
    logging.done("Folded RNA.")  # type: ignore
