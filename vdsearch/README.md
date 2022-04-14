# `vdsearch` Source Code: An Overview

`vdsearch` is built using [Typer](https://typer.tiangolo.com), a CLI framework that automatically generates CLIs from Python 3.6+ type hints.
Each command is a function in [`commands/`](/vdsearch/commands) that is imported into [`main.py`](main.py), which serves as the entry point for the CLI.
To make the CLI prettier and easier to use (_e.g._ "did you mean" when a command is misspelled), the base Typer class is overridden in [`rich_wrapper.py`](rich_wrapper.py).

Some of the commands are accelerated using [Nim](https://nim-lang.org/), a statically typed, Python-esque language.
The Nim source code is in [`nim/`](/nim).
It should be readable without too much strain by anyone familiar with Python.

To make life easier for users, some static data for RNAmotif is packaged in [`data/`](/vdsearch/data).
Any future (small scale) static data that would be useful for the pipeline should be stored here.

Finally, there are additional subcommands in [`internal/`](/vdsearch/internal) that are not exposed to the user via the CLI unless they prefix them with the `internal` command.
They mostly have to do with IO but are not part of the main pipeline.
