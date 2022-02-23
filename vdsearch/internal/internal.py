from typing import Optional
import typer
from pathlib import Path
import pandas as pd

from ..rich_wrapper import MyTyper

app = MyTyper(hidden=True)


@app.command()
def dbn2tsv(
    dbn: Path = typer.Argument(
        ...,
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="Dot bracket notation structure file",
    ),
    outfile: Optional[Path] = typer.Option(
        None,
        "--outfile",
        "-o",
        dir_okay=False,
        file_okay=True,
        help="Path to output file",
    ),
):
    """
    Convert an RNAfold .dbn file to a .tsv file.

    If no output file is specified, nothing is written to disk.
    Why? So that the function can by used in a notebook directly.
    """
    result = []

    seq_id, structure, seq, mfe = "", "", "", 0.0
    with dbn.open() as f:
        for i, line in enumerate(f):
            if i % 3 == 0:
                seq_id = line.strip()[1:]
            elif i % 3 == 1:
                seq = line.strip()
            elif i % 3 == 2:
                structure = line.strip()[0 : line.strip().rfind("(")]
                mfe = float(line.strip().split("(")[-1].split(")")[0])
                paired_percent = (len(structure) - structure.count(".")) / len(
                    structure
                )

                # count how many hairpins there are in the structure
                # a hairpin is defined by going from a ( to a )
                hairpins = 0
                last_char = "("
                for char in structure:
                    if char == "(" and last_char == "(":
                        continue
                    elif char == ")" and last_char == ")":
                        continue
                    elif char == ")" and last_char == "(":
                        hairpins += 1
                        last_char = char
                    elif char == "(" and last_char == ")":
                        last_char = char

                # the i + 1 / 3 is to account for the three rows per sequence
                # the - 1 is to make it zero indexed
                result.append(
                    pd.DataFrame(
                        [[seq_id, structure, seq, mfe, paired_percent, hairpins]],
                        index=[((i + 1) / 3) - 1],
                        columns=[
                            "seq_id",
                            "structure",
                            "seq",
                            "mfe",
                            "paired_percent",
                            "hairpins",
                        ],
                    )
                )

    # either write to disk or return the result, depending on the arguments
    result_df: pd.DataFrame = pd.concat(result)
    if outfile is not None:
        result_df.to_csv(outfile, sep="\t", index=False)
    else:
        return result_df


@app.callback()
def callback():
    """
    Unstable internal commands.

    ## Beware!

    You are entering the wilderness.
    Beyond the cozy confines of the main pipeline lies mystery, adventure, and, perhaps, danger.

    These scripts are not meant to be run by humans, only Benjamin (and maybe Neri).
    It's very likely that they are buggy, undocumented, or unrelated to the main pipeline.

    If you're not Benjamin, you should not be here.
    If you're not Neri, you should not be here.
    If you're not either, you should not be here.
    """
    pass


if __name__ == "__main__":
    app()
