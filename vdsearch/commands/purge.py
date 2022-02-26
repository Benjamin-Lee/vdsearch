import logging
from pathlib import Path
import shutil
import typer


def purge(
    directory: Path = typer.Argument(
        Path(typer.get_app_dir("vdsearch")) / "data", file_okay=False, dir_okay=True
    )
):
    """
    Purge the covariance matrix and sequence cache of all data.
    """
    data_dir = Path(typer.get_app_dir("vdsearch")) / "data"
    if data_dir.exists():
        logging.info(f"Purging data cache {data_dir}")
        shutil.rmtree(data_dir)
    else:
        logging.info("No data cache directory found.")
    logging.done("Purged data cache.")  # type: ignore
