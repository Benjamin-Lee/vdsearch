import logging
from pathlib import Path
from typing import Set

import pandas as pd
import typer


def ribozyme_filter(infernal_tsv: Path):
    ribozymes = pd.read_csv(
        infernal_tsv,
        delim_whitespace=True,
        comment="#",
        usecols=[0, 1, 2, 7, 8, 9, 15, 16],
        header=None,
        names=[
            "seq_id",
            "accession",
            "ribozyme",
            "from",
            "to",
            "strand",
            "evalue",
            "inc",
        ],
    )
    if ribozymes.shape[0] > 0:
        logging.info(
            f"Analyzing {ribozymes.shape[0]} ribozymes in {ribozymes.seq_id.unique().shape[0]} sequences to find viroid-like sequences..."
        )
    else:
        logging.done("No ribozymes present to analyze.")  # type: ignore
        return

    double_rz_ids: Set[str] = set(
        ribozymes.query("strand == '+' & evalue < 0.1").seq_id
    ) & set(ribozymes.query("strand == '-' & evalue < 0.1").seq_id)
    single_rz_ids: Set[str] = (
        set(ribozymes.query("evalue < 0.01").seq_id) - double_rz_ids
    )
    ribozy_likes_ids = double_rz_ids | single_rz_ids  # union of the two

    if len(ribozy_likes_ids) == 0:
        logging.done("No viroid-like sequences found by ribozyme search.")  # type: ignore
        return

    logging.done(  # type: ignore
        f"Found {len(ribozy_likes_ids)} viroid-like sequences. "
        f"{len(single_rz_ids)} with one ribozyme, {len(double_rz_ids)} with two ribozymes."
    )

    logging.debug("Generating output dataframes...")

    # add categorical information about how many ribozymes are in the sequence
    ribozymes.loc[ribozymes["seq_id"].isin(double_rz_ids), "Polarity"] = "(+) and (-)"
    ribozymes.loc[ribozymes["seq_id"].isin(single_rz_ids), "Polarity"] = "(+)"

    single_rzs = ribozymes.loc[ribozymes.seq_id.isin(single_rz_ids)]
    double_rzs = ribozymes.loc[ribozymes.seq_id.isin(double_rz_ids)]
    ribozy_likes = ribozymes.loc[ribozymes.seq_id.isin(ribozy_likes_ids)]

    return {
        "single_rzs": single_rzs,
        "double_rzs": double_rzs,
        "ribozy_likes": ribozy_likes,
    }


# We need to use a wrapper since, for some reason, returning values causes their return values to be printed
def ribozyme_filter_wrapper(
    infernal_tsv: Path = typer.Argument(
        ...,
        file_okay=True,
        dir_okay=False,
        exists=True,
        readable=True,
        help="Path to Infernal tabular output",
    )
):
    """Using ribozyme search results, find viroid-like sequences."""
    logging.info(f"Reading Infernal tabular output from {infernal_tsv}...")
    ribozyme_filter(infernal_tsv)
