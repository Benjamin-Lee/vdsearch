import logging
import time
from pathlib import Path
from typing import List

import typer
from rich.logging import RichHandler

from . import commands
from .internal import internal
from .rich_wrapper import MyTyper
from .types import FASTA

# configure logging for the project to use Rich
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, markup=True)],
)

# define the CLI app to use the custom Rich-enabled typer
app = MyTyper()

app.command()(commands.download_cms)
app.command()(commands.download_viroiddb)
app.command()(commands.find_circs)
app.command()(commands.dedup)
app.command()(commands.cluster)
app.command()(commands.easy_search)


# Support for using --version
# See https://github.com/tiangolo/typer/issues/52
__version__ = "0.0.1"


def version_callback(value: bool):
    if value:
        typer.echo(f"vdsearch {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Show the version and exit.",
    )
):
    """
    Workflows and utilities for finding and analyzing viroid-like circular RNAs.
    """
    return


app.add_typer(internal.app, name="internal")
if __name__ == "__main__":
    app(prog_name="vdsearch")
