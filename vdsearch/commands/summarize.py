import hashlib
import logging
from pathlib import Path

import pandas as pd
import skbio
import typer
import rich.progress
from vdsearch.types import FASTA
from vdsearch.internal import dbn2tsv
from vdsearch.utils import typer_unpacker


@typer_unpacker
def summarize(
    fasta: Path = FASTA,
    ribozyme_tsv: Path = typer.Argument(..., help="Path to ribozyme TSV file"),
    viroiddb_tsv: Path = typer.Argument(..., help="Path to ViroidDB TSV file"),
    dbn_plus: Path = typer.Argument(
        ..., help="Path to structure dbn file for + strand"
    ),
    dbn_minus: Path = typer.Argument(
        ..., help="Path to structure dbn file for - strand"
    ),
    # raw_fasta: Path = FASTA,
    source: str = typer.Argument(..., help="Source of the sequences"),
    outfile: Path = typer.Argument(..., help="Path to output file"),
    circ_tsv: Path = typer.Option(
        None, help="Path to circularity TSV file", file_okay=True, dir_okay=False
    ),
    header: bool = typer.Option(
        False, help="Include a header. Only useful for the first run"
    ),
):
    """
    Summarize the results of the analysis.
    """

    sequences = {}

    ribozymes_df = pd.read_csv(ribozyme_tsv, sep="\t")
    viroiddb_df = pd.read_csv(
        viroiddb_tsv,
        sep="\t",
        names="query,match_id,match_pident,match_alnlen,match_mismatch,match_gapopen,match_qstart,match_qend,match_tstart,match_tend,match_evalue,match_bits,match_theader,match_qcov,match_tcov,cigar".split(
            ","
        ),
    )

    # read and parse structure files
    dbn_plus_df = dbn2tsv(dbn_plus)
    dbn_plus_df.rename(
        columns={
            "mfe": "mfe_plus",
            "paired_percent": "paired_percent_plus",
            "hairpins": "hairpins_plus",
            "structure": "structure_plus",
        },
        inplace=True,
    )
    dbn_minus_df = dbn2tsv(dbn_minus)
    dbn_minus_df.rename(
        columns={
            "mfe": "mfe_minus",
            "paired_percent": "paired_percent_minus",
            "hairpins": "hairpins_minus",
            "structure": "structure_minus",
        },
        inplace=True,
    )

    if circ_tsv is not None and circ_tsv.exists():
        circ_df = pd.read_csv(circ_tsv, sep="\t")

    with rich.progress.open(fasta, transient=True) as f:  # type: ignore
        for seq in skbio.read(f, format="fasta", constructor=skbio.DNA):
            seq_data = {}
            seq_id = seq.metadata["id"]
            # basic sequence information
            seq_data["seq_id"] = seq_id
            seq_data["description"] = seq.metadata["description"]
            seq_data["unit_length"] = len(seq)
            seq_data["gc_content"] = seq.gc_content()

            # ribozyme information
            seq_ribozymes = ribozymes_df.query(f"seq_id == '{seq_id}'").sort_values(
                by=["evalue"], ascending=True
            )
            seq_data["has_ribozymes"] = (
                bool(len(seq_ribozymes))
                and ["Pospi_RY"] != seq_ribozymes.ribozyme.unique().tolist()
            )

            # I've renamed Polarity to symmetry (and made it boolean) but I'm keeping backwards compatibility
            if not seq_data["has_ribozymes"]:
                seq_data["symmetric"] = False
            elif "Polarity" in seq_ribozymes.columns:
                seq_data["symmetric"] = (
                    seq_ribozymes.Polarity.unique()[0] == "(+) and (-)"
                )
            else:
                seq_data["symmetric"] = seq_ribozymes.symmetric.unique()[0]

            if seq_data["has_ribozymes"]:
                rzs_present = [
                    seq_ribozymes.loc[seq_ribozymes.strand == strand]
                    for strand in ["+", "-"]
                ]
                use_plus_rz = True
                use_minus_rz = True

                # if the sequence isn't symmetric, we only want to look at the strand with the best ribozyme
                # note that we want to use evalue not bitscore
                if not seq_data["symmetric"]:
                    best_plus = rzs_present[0].evalue.min()
                    best_minus = rzs_present[1].evalue.min()
                    if pd.isna(best_minus) or best_plus < best_minus:
                        use_plus_rz = True
                        use_minus_rz = False
                    else:
                        use_plus_rz = False
                        use_minus_rz = True

                if rzs_present[0].shape[0] and use_plus_rz:
                    seq_data["rz_plus"] = rzs_present[0].ribozyme.values[0]
                    seq_data["rz_plus_evalue"] = rzs_present[0].evalue.values[0]
                    seq_data["rz_plus_score"] = rzs_present[0].score.values[0]
                    seq_data["rz_plus_from"] = rzs_present[0]["from"].values[0]
                    seq_data["rz_plus_to"] = rzs_present[0]["to"].values[0]
                if rzs_present[1].shape[0] and use_minus_rz:
                    seq_data["rz_minus"] = rzs_present[1].ribozyme.values[0]
                    seq_data["rz_minus_evalue"] = rzs_present[1].evalue.values[0]
                    seq_data["rz_minus_score"] = rzs_present[1].score.values[0]
                    seq_data["rz_minus_from"] = rzs_present[1]["from"].values[0]
                    seq_data["rz_minus_to"] = rzs_present[1]["to"].values[0]

            # viroiddb information
            viroiddb_match = (
                viroiddb_df.query(f"query == '{seq_id}'")
                .sort_values(by="match_bits", ascending=False)
                .drop_duplicates(subset="query", keep="first")
                .to_dict(orient="records")
            )
            if len(viroiddb_match):
                seq_data.update((viroiddb_match[0]))

            # merge in the structure data
            seq_data.update(
                dbn_plus_df.query(f"seq_id == '{seq_id}'").to_dict(orient="records")[0]
            )
            seq_data.update(
                dbn_minus_df.query(f"seq_id == '{seq_id}'").to_dict(orient="records")[0]
            )

            # merge in the circ data
            if circ_tsv.exists() and circ_df is not None:
                seq_data.update(
                    circ_df.query(f"seq_id == '{seq_id}'").to_dict(orient="records")[0]
                )
            else:
                seq_data["original_length"] = None
                seq_data["ratio"] = None

            seq_data["source"] = source

            # generate a new id based on the md5 of the canonicalized sequence
            new_id = (
                "NV_"
                + hashlib.blake2b(str(seq).encode("utf-8"), digest_size=8).hexdigest()
            )
            seq_data["vdsearch_id"] = new_id

            # add the sequence to the output
            sequences[new_id] = seq_data

    res = pd.DataFrame.from_dict(sequences, orient="index")[
        [
            "vdsearch_id",
            "seq_id",
            "description",
            "source",
            "unit_length",
            "original_length",
            "ratio",
            "gc_content",
            # ribozyme info
            "has_ribozymes",
            "symmetric",
            "rz_plus",
            "rz_minus",
            "rz_plus_from",
            "rz_plus_to",
            "rz_plus_evalue",
            "rz_plus_score",
            "rz_minus_from",
            "rz_minus_to",
            "rz_minus_evalue",
            "rz_minus_score",
            # structure info
            "mfe_plus",
            "paired_percent_plus",
            "hairpins_plus",
            "mfe_minus",
            "paired_percent_minus",
            "hairpins_minus",
            # vdsearch results
            "match_id",
            "match_pident",
            "match_alnlen",
            "match_mismatch",
            "match_gapopen",
            "match_qstart",
            "match_qend",
            "match_tstart",
            "match_tend",
            "match_evalue",
            "match_bits",
            "match_theader",
            "match_qcov",
            "match_tcov",
            "seq",
            "structure_plus",
            "structure_minus",
        ]
    ]
    res.to_csv(outfile, sep="\t", index=False, header=header)
    logging.done(  # type: ignore
        f"Found {len(res.vdsearch_id.to_list())} total viroid-like sequences. "
    )
