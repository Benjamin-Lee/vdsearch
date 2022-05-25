import logging
from pathlib import Path
import subprocess
from typing import Optional

import typer

from vdsearch.types import FASTA, Threads
from vdsearch.utils import check_executable_exists, typer_unpacker


@typer_unpacker
def fold(
    fasta: Path = FASTA,
    output: Optional[Path] = typer.Option(None),
    threads: int = Threads,
    ps: bool = typer.Option(False, help="Include PostScript output"),
    ps_dir: Path = typer.Option(
        None, help="Directory to store PostScript output. Implies --ps."
    ),
    temp: int = typer.Option(25, help="Number of degrees C to use for folding"),
):
    """Predict secondary structures of circRNAs.

    ## References

    > Lorenz, R., Bernhart, S.H., Höner zu Siederdissen, C., Tafer, H., Flamm, C., Stadler, P.F., et al., 2011.
    > ViennaRNA Package 2.0.
    > Algorithms Mol Biol 6, 26.
    > https://doi.org/10.1186/1748-7188-6-26

    """
    check_executable_exists("RNAfold", "ViennaRNA")

    if output is None:
        logging.debug("Using default output filename.")
        output = Path(fasta).with_suffix(".dbn")

    # we might need to make a directory for the PS files
    if ps_dir is not None and not ps_dir.exists():
        ps_dir.mkdir(parents=True)

    command = (
        f"{'cd ' + str(ps_dir) + ' &&' if ps_dir is not None else '' } "  # support for PS output to a directory
        # We need the ps or  ps_dir because ps_dir implies ps
        f"RNAfold --circ {'' if ps or ps_dir is not None else '--noPS'} --jobs={threads} --temp={temp} {fasta.absolute()} > {output.absolute()} 2> /dev/null"
    )
    logging.debug(f"{command=}")
    logging.info(f"Folding {fasta}")
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(
            f"Command [cyan]{command}[/] failed with exit code [red]{e.returncode}[/]"
        )
        raise e
    if ps_dir:
        logging.info(f"PostScript output written to {ps_dir}")
    logging.done("Folded RNA.")  # type: ignore
