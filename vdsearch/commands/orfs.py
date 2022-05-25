from pathlib import Path
from typing import List
import rich.progress
import skbio

import typer

from vdsearch.types import FASTA
from vdsearch.utils import typer_unpacker


@typer_unpacker
def orfs(
    fasta: Path = FASTA,
    outfile: Path = typer.Argument(..., help="Path to output file", dir_okay=False),
    orf_len: int = typer.Option(100, help="Minimum ORF length in amino acids"),
):
    """
    Extract the ORFs from a FASTA file containing circRNAs.
    """
    import orfipy_core

    with rich.progress.open(  # type: ignore
        fasta,
        transient=True,
        description="Extracting ORFs",
    ) as f, outfile.open("w") as o:
        for seq in skbio.read(f, format="fasta", constructor=skbio.DNA):
            # we can't deal with degenerate sequences
            if seq.has_degenerates():
                continue

            # perform the self-concatenation
            doubled = str(seq) + str(seq)

            # we use a set to deduplicate ORFs sequences resulting from self-concatenation
            orfs = set()
            unique_orfs: List[
                skbio.Protein
            ] = []  # note: unique within the sequence's ORFs

            for start, stop, strand, _ in orfipy_core.orfs(
                doubled,
                minlen=orf_len * 3,
                name=seq.metadata["id"],
                starts=["ATG"],
                stops=["TAA", "TAG", "TGA"],
            ):
                # get the translated ORF
                orf_seq = skbio.DNA(doubled[start:stop])
                if strand == "-":
                    orf_seq = orf_seq.reverse_complement()
                translated = orf_seq.translate()

                # If the ORF is unique within the sequence, add it to the sequence's ORF list
                # TODO: keep only the first ORF
                if str(translated) not in orfs:
                    orfs.add(str(translated))
                    unique_orfs.append(
                        skbio.Protein(
                            translated,
                            metadata={
                                "id": seq.metadata["id"]
                                + "_start_"
                                + str(start if strand == "+" else stop)
                                + "_stop_"
                                + str(stop if strand == "+" else start)
                            },
                        )
                    )

            # write the unique ORFs out
            for orf in unique_orfs:
                orf.write(o, format="fasta")
