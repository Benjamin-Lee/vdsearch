import logging
import time
from pathlib import Path
from typing import List

import typer


from . import commands
from .internal import internal
from .rich_wrapper import MyTyper
from .types import FASTA


# define the CLI app to use the custom Rich-enabled typer
app = MyTyper()

app.command()(commands.download_cms)  # type: ignore
app.command()(commands.download_viroiddb)  # type: ignore
app.command()(commands.find_circs)  # type: ignore
app.command()(commands.dedup)  # type: ignore
app.command()(commands.cluster)  # type: ignore
app.command()(commands.easy_search)  # type: ignore
app.command()(commands.canonicalize)  # type: ignore
app.command()(commands.infernal)  # type: ignore
app.command()(commands.ribozyme_filter)  # type: ignore


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
    ),
    verbose: bool = typer.Option(False, "--verbose", help="Show verbose output."),
):
    """
    Workflows and utilities for finding and analyzing viroid-like circular RNAs.
    """
    if verbose:
        logging.basicConfig(level=logging.DEBUG)


app.add_typer(internal.app, name="internal")
if __name__ == "__main__":
    app(prog_name="vdsearch")
