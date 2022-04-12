import logging
from pathlib import Path
import typer
import fileinput


from vdsearch.types import Threads

"""
A rudimentary URL downloader (like wget or curl) to demonstrate Rich progress bars.

Taken shamelessly from https://github.com/Textualize/rich/blob/4b3b6531ad349312f4df6e284bc1831dcf59ada6/examples/downloader.py
"""

import os.path
import sys
from concurrent.futures import as_completed, ThreadPoolExecutor
import signal
from functools import partial
from threading import Event
from typing import Iterable
from urllib.request import urlopen

from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    TaskID,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)

progress = Progress(
    TextColumn("[bold blue]{task.fields[filename]}", justify="right"),
    BarColumn(bar_width=None),
    "[progress.percentage]{task.percentage:>3.1f}%",
    "•",
    DownloadColumn(),
    "•",
    TransferSpeedColumn(),
    "•",
    TimeRemainingColumn(),
)


done_event = Event()


def handle_sigint(signum, frame):
    done_event.set()


signal.signal(signal.SIGINT, handle_sigint)


def copy_url(task_id: TaskID, url: str, path: str) -> None:
    """Copy data from a url to a local file."""
    # progress.console.log(f"Requesting {url}")
    response = urlopen(url)
    # This will break if the response doesn't contain content length
    progress.update(task_id, total=int(response.info()["Content-length"]))
    with open(path, "wb") as dest_file:
        progress.start_task(task_id)
        for data in iter(partial(response.read, 32768), b""):
            dest_file.write(data)
            progress.update(task_id, advance=len(data))
            if done_event.is_set():
                return


def download(urls: Iterable[str], dest_dir: str):
    """Download multuple files to the given directory."""

    with progress:
        with ThreadPoolExecutor(max_workers=8) as pool:
            for url in urls:
                filename = url.split("/")[-2]
                dest_path = os.path.join(dest_dir, filename)
                task_id = progress.add_task("download", filename=filename, start=False)
                pool.submit(copy_url, task_id, url, dest_path)


def download_viroiddb():
    """Download the latest ViroidDB dataset."""
    print("downloading vdb")
    # with Halo(text="Downloading ViroidDB dataset...") as sp:
    #     time.sleep(2)
    #     sp.succeed("Downloaded ViroidDB")


cms = {
    "RF03160": {"description": "type-P1 twister ribozyme", "name": "twister-P1"},
    "RF02684": {"description": "Type-P5 twister ribozyme", "name": "Twister-P5"},
    "RF02678": {"description": "Hatchet ribozyme", "name": "Hatchet"},
    "RF02681": {
        "description": "Twister sister ribozyme ",
        "name": "Twister-sister",
    },
    "RF00008": {
        "description": "Hammerhead ribozyme (type III)",
        "name": "Hammerhead_3",
    },
    "RF04190": {
        "description": "Hairpin ribozyme 1 from viruses-like metatranscriptomes",
        "name": "Hairpin-meta1",
    },
    "RF04188": {
        "description": "Hovlinc ribozyme (hominin vlincRNA-located)",
        "name": "Hovlinc",
    },
    "RF04191": {
        "description": "Hairpin ribozyme 2 from viruses-like metatranscriptomes",
        "name": "Hairpin-meta2",
    },
    "RF02682": {
        "description": "HDV ribozyme from F. prausnitzii",
        "name": "HDV-F-prausnitzii",
    },
    "RF02679": {"description": "Pistol ribozyme", "name": "Pistol"},
    "RF00094": {
        "description": "Hepatitis delta virus ribozyme",
        "name": "HDV_ribozyme",
    },
    "RF02275": {"description": "Hammerhead ribozyme HH9", "name": "Hammerhead_HH9"},
    "RF00163": {
        "description": "Hammerhead ribozyme (type I)",
        "name": "Hammerhead_1",
    },
    "RF02276": {
        "description": "Hammerhead ribozyme (type II)",
        "name": "Hammerhead_II",
    },
    "RF02277": {
        "description": "Hammerhead ribozyme HH10",
        "name": "Hammerhead_HH10",
    },
    "RF03152": {
        "description": "RAGATH-1 hammerhead ribozyme",
        "name": "RAGATH-1-hammerhead",
    },
    "RF00622": {
        "description": "Mammalian CPEB3 ribozyme",
        "name": "CPEB3_ribozyme",
    },
    "RF03154": {"description": "type-P3 twister ribozyme", "name": "twister-P3"},
    "RF01787": {"description": "drz-agam-1 ribozyme", "name": "drz-agam-1"},
    "RF00173": {"description": "Hairpin ribozyme", "name": "Hairpin"},
    "RF01788": {"description": "drz-agam-2-2 ribozyme", "name": "drz-agam-2-2"},
    "RF00362": {
        "description": "Pospiviroid RY motif stem loop",
        "name": "Pospi_RY",
    },
}


def download_cms(
    force: bool = typer.Option(False, help="Force download of the reference CMs"),
    threads: int = Threads,
):
    """Download the latest covariance matrices for self-cleaving ribozymes."""

    app_dir = Path(typer.get_app_dir("vdsearch"))
    cms_dir = app_dir / "data" / "cms" / "raw"
    cms_dir.mkdir(exist_ok=True, parents=True)

    needed_cms = [k for k in cms.keys() if not (cms_dir / f"{k}").exists()]
    if needed_cms or force:
        logging.info("Downloading ribozyme CMs from Rfam...")
        download(
            [
                f"https://rfam.org/family/{k}/cm"
                for k in (needed_cms if not force else cms.keys())
            ],
            str(cms_dir),
        )
        logging.done("Downloaded ribozyme CMs")  # type: ignore
    else:
        logging.done("All ribozyme CMs are already downloaded.")  # type: ignore

    logging.debug("Verifying CMs...")
    for k, v in cms.items():
        cms_file = cms_dir / f"{k}"

        # check if the file exists
        if not cms_file.exists():
            raise FileNotFoundError(f"Could not find {cms_file}")

        lines = cms_file.read_text().splitlines()

        # check that the name is correct
        name = lines[1].split()[1]
        if name != v["name"]:
            raise ValueError(f"{cms_file} has name {name} instead of {v['name']}")

        # check that the file is complete
        if not lines.count("//") == 2:
            raise ValueError(f"{cms_file.as_posix()} does not contain two // lines")
    logging.debug("CMs verified")

    logging.debug("Merging CMs into one file...")
    merged = app_dir / "data" / "cms" / "merged"
    merged.parent.mkdir(exist_ok=True, parents=True)
    with open(merged, "w") as fout, fileinput.input(
        [cms_dir / cm for cm in cms.keys()]
    ) as fin:
        for line in fin:
            fout.write(line)
    logging.debug("Merged CMs")

    logging.debug("Verifying merged CMs...")
    if merged.read_text().splitlines().count("//") != len(cms) * 2:
        raise ValueError(f"{merged.as_posix()} does not contain 2*{len(cms)} // lines")
    logging.debug("Verified merged CMs")
