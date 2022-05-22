from enum import Enum
import logging
from pathlib import Path
from typing import Dict, Literal, Optional

import pandas as pd
import typer

from vdsearch.types import ReferenceCms


class CmCutoffType(str, Enum):
    GA = "GA"
    NC = "NC"
    TC = "TC"


def parse_cm_file(path: Path) -> Dict[str, Dict[str, float]]:
    """Given a path to a CM file, parse it and return a dictionary of the form:

    ```
    {
        'Twister-P5': {
            'GA': 45.0,
            'TC': 45.3,
            'NC': 32.3,
    }
    `
    """
    cutoffs = {}

    with path.open() as f:
        last_name = ""
        last_cutoff = {
            "GA": 0.0,
            "TC": 0.0,
            "NC": 0.0,
        }

        for line in f:
            if "INFERNAL" in line:
                if last_name:
                    cutoffs[last_name] = last_cutoff

                last_name = ""
                last_cutoff = {
                    "GA": 0.0,
                    "TC": 0.0,
                    "NC": 0.0,
                }
            elif line.startswith("NAME"):
                last_name = line.split()[1]
            elif line.startswith("GA"):
                last_cutoff["GA"] = float(line.split()[1])
            elif line.startswith("TC"):
                last_cutoff["TC"] = float(line.split()[1])
            elif line.startswith("NC"):
                last_cutoff["NC"] = float(line.split()[1])

        # at the end, we need to add the last one
        cutoffs[last_name] = last_cutoff
    return cutoffs


def ribozyme_filter(
    infernal_tblout: Path,
    rnamotif_txt: Path = None,
    rnamotif_name: str = None,
    output_tsv: Optional[Path] = None,
    use_cm_cutoff: bool = True,
    cm_file: Optional[Path] = None,
    cm_cutoff_type: Literal["GA", "TC", "NC"] = "GA",
    use_evalue_cutoff: bool = True,
    max_evalue: float = 0.01,
):
    ribozymes = pd.read_csv(
        infernal_tblout,
        delim_whitespace=True,
        comment="#",
        usecols=[0, 1, 2, 7, 8, 9, 14, 15, 16],
        header=None,
        names=[
            "seq_id",
            "accession",
            "ribozyme",
            "from",
            "to",
            "strand",
            "score",
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

    cutoffs = parse_cm_file(cm_file) if use_cm_cutoff and cm_file else {}

    rz_plus = set()
    rz_minus = set()
    rz_significant = set()

    for rz_name, rz_df in ribozymes.groupby(["ribozyme"]):

        # short circuit if we aren't doing verbose output
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            logging.debug(
                f"Analyzing {rz_name=}. "
                f"{rz_df.shape[0]} sequences, "
                f"{rz_df.query(f'evalue < {max_evalue}').shape[0]} significant."
            )

        # all seqs with ribozymes above cutoffs are counted as significant
        if use_cm_cutoff and cutoffs:
            rz_plus.update(
                rz_df.query(
                    f"strand == '+' & score > {cutoffs[rz_name][cm_cutoff_type]}"
                ).seq_id
            )
            rz_minus.update(
                rz_df.query(
                    f"strand == '-' & score > {cutoffs[rz_name][cm_cutoff_type]}"
                ).seq_id
            )
            # signficant ribozymes are either above cutoff with extra low evalue
            rz_significant.update(
                rz_df.query(f"score > {cutoffs[rz_name][cm_cutoff_type]}").seq_id
            )

        if use_evalue_cutoff:
            # we will also use the evalue cutoff to determine if a ribozyme is present
            rz_plus.update(
                rz_df.query(f"strand == '+' & evalue < {max_evalue ** 0.5}").seq_id
            )
            rz_minus.update(
                rz_df.query(f"strand == '-' & evalue < {max_evalue ** 0.5}").seq_id
            )
            rz_significant.update(rz_df.query(f"evalue < {max_evalue}").seq_id)

    # Parse and add RNAmotif hits
    if rnamotif_txt and rnamotif_txt.exists():
        rnamotifs = pd.read_csv(
            rnamotif_txt,
            delim_whitespace=True,
            comment="#",
            header=None,
            names=["seq_id", "score", "strand", "from_", "length"],
            usecols=[0, 1, 2, 3, 4],
        )
        rz_plus.update(rnamotifs.query("strand == 0").seq_id)
        rz_minus.update(rnamotifs.query("strand == 1").seq_id)
        # note that there's no add to rz_significant here
        # this is since I'm not sure if the RNAmotif hits are significant

        # merge the rnamotifs df into the ribozymes df
        # iterate over rnamotifs
        for rnamotif in rnamotifs.itertuples(index=False):
            ribozymes = pd.concat(
                [
                    ribozymes,
                    pd.DataFrame(
                        [
                            {
                                "seq_id": rnamotif.seq_id,
                                "ribozyme": rnamotif_name,
                                "strand": "+" if rnamotif.strand == 0 else "-",
                                "evalue": max_evalue**0.5,
                                "from": rnamotif.from_,
                                "to": rnamotif.from_
                                + (
                                    rnamotif.length * 1
                                    if rnamotif.strand == 0
                                    else -1 * rnamotif.length
                                ),
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )

    # any sequence with a ribozyme match (even weaker than cutoff) is counted as significant if there are two
    double_rz_ids = rz_plus.intersection(rz_minus)
    # all sequences with a significant ribozyme
    single_rz_ids = rz_significant.difference(double_rz_ids)
    # union of the two
    ribozy_likes_ids = double_rz_ids | single_rz_ids

    if len(ribozy_likes_ids) == 0:
        logging.done("No viroid-like sequences found by ribozyme search.")  # type: ignore
        return

    logging.done(  # type: ignore
        f"Found {len(ribozy_likes_ids)} viroid-like sequences with ribozymes. "
        f"{len(single_rz_ids)} with one ribozyme, {len(double_rz_ids)} with two ribozymes."
    )

    logging.debug("Generating output dataframes...")

    # add categorical information about how many ribozymes are in the sequence
    ribozymes.loc[ribozymes["seq_id"].isin(double_rz_ids), "symmetric"] = True
    ribozymes.loc[ribozymes["seq_id"].isin(single_rz_ids), "symmetric"] = False

    single_rzs = ribozymes.loc[ribozymes.seq_id.isin(single_rz_ids)]
    double_rzs = ribozymes.loc[ribozymes.seq_id.isin(double_rz_ids)]
    ribozy_likes = ribozymes.loc[ribozymes.seq_id.isin(ribozy_likes_ids)]

    if ribozy_likes.shape[0] and output_tsv:
        ribozy_likes.sort_values(by="evalue").to_csv(
            output_tsv,
            sep="\t",
            index=False,
        )

    return {
        "single_rzs": single_rzs,
        "double_rzs": double_rzs,
        "ribozy_likes": ribozy_likes,
    }


# We need to use a wrapper since, for some reason, returning values causes their return values to be printed
def ribozyme_filter_wrapper(
    infernal_tblout: Path = typer.Argument(
        ...,
        file_okay=True,
        dir_okay=False,
        exists=True,
        readable=True,
        help="Path to Infernal tabular output",
    ),
    rnamotif_txt: Path = typer.Option(
        None,
        help="Path to pruned RNAmotif .txt file for viroid-like sequences. Note that e-values will be fixed at the square root of --max-evalue. Concatenated outputs of RNAmotif are not supported.",
    ),
    rnamotif_name: str = typer.Option(
        None, help="Name of motif to use in the ribozyme column of the output"
    ),
    output_tsv: Path = typer.Option(
        None,
        help="Path to output TSV file with Infernal results for viroid-like sequences",
    ),
    use_cm_cutoff: bool = typer.Option(
        True, help="Use CM cutoffs to determine if a ribozyme is present."
    ),
    cm_file: Optional[Path] = ReferenceCms,
    cm_cutoff_type: CmCutoffType = typer.Option("GA"),
    use_evalue_cutoff: bool = typer.Option(
        True, help="Use evalue cutoff to determine if a ribozyme is present."
    ),
    max_evalue: float = typer.Option(
        0.01, help="Maximum evalue to use when determining if a ribozyme is present."
    ),
):
    """Using ribozyme search results, find viroid-like sequences."""
    logging.info(f"Reading Infernal tabular output from {infernal_tblout}")
    ribozyme_filter(
        infernal_tblout,
        rnamotif_txt=rnamotif_txt,
        rnamotif_name=rnamotif_name,
        output_tsv=output_tsv,
        use_cm_cutoff=use_cm_cutoff,
        cm_file=cm_file,
        cm_cutoff_type=cm_cutoff_type.value,  # type: ignore
        use_evalue_cutoff=use_evalue_cutoff,
        max_evalue=max_evalue,
    )
    logging.done(f"Wrote ribozyme data for viroid-like sequences to {output_tsv}")  # type: ignore

    # table = rich.table.Table(
    #     highlight=True,
    #     title="Ribozyme Search Results",
    #     box=rich.box.ROUNDED,
    #     show_footer=True,
    # )
    # table.add_column("Ribozyme", style="magenta", footer="Total")
    # table.add_column("(+) only count", footer=f"{len(results['single_rzs'])}")
    # table.add_column("(+) and (-) count", footer=f"{len(results['double_rzs'])}")

    # for ribozyme, ribozyme_df in results["ribozy_likes"].groupby(["ribozyme"]):
    #     table.add_row(
    #         f"{ribozyme}",
    #         str(ribozyme_df.query("Polarity == '(+)'").seq_id.nunique()),
    #         str(ribozyme_df.query("Polarity == '(+) and (-)'").seq_id.nunique()),
    #     )
    # rich.get_console().log(table)
