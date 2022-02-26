import functools
import logging
import shutil
from typing import Callable

import click
from typer.models import ParameterInfo


def check_executable_exists(
    name: str, display_name: str = None, test_command: str = None
):
    """
    Check if an executable exists in the system path.
    """
    if shutil.which(name) is None:
        logging.error(f"Unable to find {display_name} executable.")
        raise click.Abort(
            f"Unable to find {display_name or name} executable. Are you sure that it's installed? "
            f"To test if it is, run the following command: {test_command or name}"
        )


def typer_unpacker(f: Callable):
    """
    Make a Typer function into a normal function.

     from https://github.com/tiangolo/typer/issues/279#issuecomment-841875218
    """

    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        default_values = f.__defaults__
        patched_defaults = [
            value.default if isinstance(value, ParameterInfo) else value
            for value in default_values
        ]
        f.__defaults__ = tuple(patched_defaults)

        return f(*args, **kwargs)

    return wrapper
