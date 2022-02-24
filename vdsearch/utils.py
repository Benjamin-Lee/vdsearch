import logging
import shutil

import click


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
