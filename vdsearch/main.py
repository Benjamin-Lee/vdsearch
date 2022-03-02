import logging

import rich
import typer
from rich.logging import RichHandler

from . import commands
from .internal import internal
from .rich_wrapper import MyTyper


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
app.command("ribozyme-filter")(commands.ribozyme_filter_wrapper)  # type: ignore
app.command()(commands.purge)  # type: ignore
app.command()(commands.write_seqs)  # type: ignore

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
    level = logging.DEBUG if verbose else logging.INFO
    # configure logging for the project to use Rich
    console = rich.console.Console(
        theme=rich.theme.Theme(
            {"logging.level.done": "green", "logging.level.debug": "dim"}
        )
    )
    logging.basicConfig(
        level=level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True, markup=True, console=console)],
    )


app.add_typer(internal.app, name="internal")
if __name__ == "__main__":
    app(prog_name="vdsearch")
