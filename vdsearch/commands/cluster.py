import itertools
import logging
import random
import shutil
import subprocess
from collections import defaultdict
from enum import Enum
from pathlib import Path
import time

import click
import igraph as ig
import numpy as np
import pandas as pd
import rich.progress
import typer
from vdsearch.types import FASTA, Threads
from vdsearch.utils import check_executable_exists, typer_unpacker


class PRESET(str, Enum):
    DEFAULT = "none"
    NT_PRECLUSTER = "nt-precluster"
    NT_CLUSTER = "nt-cluster"
    ORF = "orfs"


NT_CLUSTER_COLNAMES = (
    "query,target,pident,qlen,tlen,qstart,qend,tstart,tend,alnlen,evalue,bits,qcov,tcov"
)


def pick_resolution(graph, target_avg_ani, steps=101, seed=1953):
    """
    Given a graph (`graph`) and a target average within-cluster ANI (`target_avg_ani`),
    get the resolution parameter that will provide the closest within-community average
    edge weight using the Leiden clustering algorithm.
    """
    # The `last_res` variable will store the resolution (`res`) value of the previous
    # iteration, while the `last_avg_weight` will store the difference between the average
    # AAI of the previous iteration and the target average ANI (`target_avg_ani`)
    last_res = None
    last_avg_weight = None
    # Iterate through resolution parameters (`res`)
    for res in rich.progress.track(
        np.linspace(1, 0, steps),
        description="Picking resolution",
        total=steps,
        transient=True,
    ):
        random.seed(seed)
        # Find the communities using the current resolution parameter
        communities = graph.community_leiden(weights="weight", resolution_parameter=res)
        # Initiate the list that will store the edges' ANI
        ani_list = []
        # Iterate through the communities subgraphs ordered by length (largest → smallest)
        for sg in reversed(np.argsort(communities.sizes())):
            members = [i.attributes()["name"] for i in communities.subgraph(sg).vs]
            # If the community has a single member, break the loop
            if len(members) < 2:
                break
            # Add all the pairwise ANI value to the `ani_list` list
            for i, j in itertools.combinations(members, 2):
                try:
                    edge = graph.get_eid(i, j)
                    edge = graph.es[edge]
                    ani_list.append(edge.attributes()["ani"])
                # In case two nodes are not connected, skip the pair
                except Exception:
                    continue
        # Compute the average value of all the edges in the iteration
        current_avg_ani = np.mean(ani_list) if len(ani_list) else 1
        # If this is the first iteration (that is, `res == 1`)
        if res == 1:
            last_res = res
            last_avg_weight = current_avg_ani
        # If the difference between the curent iteration average ANI and the target ANI
        # is less than the difference of the previous iteration, store the value and keep
        # iterating through `res`
        elif np.abs(target_avg_ani - current_avg_ani) <= np.abs(
            target_avg_ani - last_avg_weight
        ):
            last_res = res
            last_avg_weight = current_avg_ani
        # If the difference between the curent iteration average ANI and the target ANI
        # is greater than the last iteration, break the loop
        else:
            break
    return last_res


@typer_unpacker
def AvA2cluster(
    path: Path = typer.Argument(..., help="Path to the MMseqs all-vs-all output file"),
    outfile: Path = typer.Argument(..., help="Path to the output file"),
    ani: float = typer.Option(0.9, help="Target average ANI", min=0, max=1),
    min_cov: float = typer.Option(0.5, help="Minimum coverage.", min=0, max=1),
    columns: str = typer.Option(NT_CLUSTER_COLNAMES, help="Columns in the input file"),
):
    """Convert an MMseqs all-vs-all search output file to a cluster file.

    ## Notes

    This method assumes that the all-vs-all search was conducted on sequences concatenated to themselves to account for circularity.
    It won't work correctly if the input sequences are linear.

    ## References

    This method is based on the following research:

    1. Traag, V. A., Waltman. L., van Eck, N. J., 2019.
       From Louvain to Leiden: guaranteeing well-connected communities.
       Scientific Reports 9, 5233.
       <https://doi.org/10.1038/s41598-019-41695-z>
    2. Csárdi, G., Nepusz, T., 2006.
       The igraph software package for complex network research.
       InterJournal Complex Systems, 1695.
    """
    # we need special logic for performing ANI clustering
    logging.info(
        f"Converting all-vs-all search to clusters with target average {ani*100:.0f}% ANI and {min_cov*100:.0f}% minimum coverage..."
    )
    logging.debug("Reading all-vs-all search file")
    df = pd.read_csv(path, sep="\t", names=columns.split(","))

    # Cull each {query,target} pair to their best available alignment
    logging.debug("Sorting all-vs-all search file by bit score")
    df.sort_values(by=["query", "target", "bits"], inplace=True, ascending=False)
    logging.debug("Dropping duplicate alignments, keeping the best")
    df.drop_duplicates(subset=["query", "target"], keep="first", inplace=True)

    # for each pair, the length of the shorter sequence and divide by two to deconcatenate
    shorter = np.minimum(df.qlen, df.tlen) / 2

    # we cap the aligned length to the shorter of the sequences
    alnlen = np.minimum(shorter, df.alnlen)

    df["ANI"] = (df.pident * alnlen) / (shorter)
    df["ANI"] /= 100

    # since the sequences were duplicated, the coverages must be doubled
    #
    # if the covered region was more than 0.5, we have a longer-than-unit
    # alignment so we have to cap it at 1.0
    df["AF"] = np.minimum(np.minimum(df.qcov * 2, df.tcov * 2), 1.0)

    # Initiate the lists that will store the edges between two sequences and
    # their respective weights (ANI × coverage) and ANI
    nodes = []
    edges = []
    weights = defaultdict(list)
    anis = defaultdict(list)

    for row in rich.progress.track(
        df.itertuples(), total=len(df), transient=True, description="Generating edges"
    ):
        nodes.extend([row.query, row.target])
        if row.AF > min_cov:
            pair = tuple(sorted([row.query, row.target]))
            edges.append(pair)
            anis[pair].append(row.ANI)
            weight = row.ANI  # * max([row.qcov, row.tcov])
            weights[pair].append(weight)

    # Remove duplicated nodes
    nodes = sorted(set(nodes))
    # Take the mean ANI and weight for each pair of sequences
    anis = [np.mean(anis[pair]) for pair in edges]
    weights = [np.mean(weights[pair]) for pair in edges]
    # Create a graph from the node list and add edges weighted by the ANI between
    # the sequences being connected
    graph = ig.Graph()
    logging.debug("Adding nodes to graph")
    graph.add_vertices(nodes)
    logging.debug("Adding edges to graph")
    graph.add_edges(edges, attributes={"weight": weights, "ani": weights})
    logging.info("Picking the best resolution for Leiden clustering")
    leiden_resolution = pick_resolution(graph, target_avg_ani=ani)
    clusters = graph.community_leiden(
        weights="weight", resolution_parameter=leiden_resolution
    )
    with open(outfile, "w") as fout:
        for i in reversed(np.argsort(clusters.sizes())):
            subgraph = clusters.subgraph(i)
            members = [v.attributes()["name"] for v in subgraph.vs]
            clu_rep = members[0]  # assign the representative to the first member
            for member in members:
                fout.write(f"{clu_rep}\t{member}\n")
    logging.done("Done postprocessing.")  # type: ignore


@typer_unpacker
def cluster(
    fasta: Path = FASTA,
    prefix: str = typer.Option(None, help="Prefix for the output files"),
    tmpdir: Path = typer.Option(
        Path(f"tmp.{int(time.time())}"),
        help="Path to temporary directory to use for intermediate files",
    ),
    preset: PRESET = typer.Option(
        "none",
        help="Mode to use for clustering. Options are 'nt-precluster' or 'orfs'",
    ),
    lin: bool = typer.Option(
        False,
        help="Use MMseqs linear commands (*i.e.* `easy-linclust` and `easy-linsearch`). May be faster but less accurate.",
    ),
    threads: int = Threads,
    min_seq_id=typer.Option(None, hidden=True),
    min_aln_len=typer.Option(None, hidden=True),
    seq_id_mode=typer.Option(None, hidden=True),
    cov_mode=typer.Option(None, hidden=True),
    coverage=typer.Option(None, hidden=True),
    evalue=typer.Option(None, hidden=True),
    sensitivity=typer.Option(None, hidden=True),
    cluster_mode=typer.Option(None, hidden=True),
    max_seqs=typer.Option(None, hidden=True),
    k=typer.Option(None, hidden=True),
):
    """Cluster circRNAs or ORFs.

    This command will cluster the sequences in the input FASTA file using one of several methods.
    The default is to just call `MMseqs2 easy-cluster` with its default parameters.
    There are also multiple presets that can be used specifically to cluster circRNAs or ORFs.

    ## Presets

    The `--preset` option can be used to specify a preset mode. The following presets are available:
    - `nt-precluster`: Cluster sequences using very strict settings. Useful for removing near-identical sequences.
    - `nt-cluster`: Compute an all-vs-all distance matrix for the circRNAs in a circular-aware fashion and then use the Leiden algorithm to produce ANI90 clusters.
      When this preset is used, the `ava2cluster` command will be called to postprocess the distance matrix into a MMseqs2 cluster TSV file.
    - `orfs`: Cluster ORFs.

    ## References

    This method is based on the following research:

    1. Steinegger, M., Soeding, J., 2017.
       MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.
       Nature Biotechnology 35, pages 1026–1028
       <https://doi.org/10.1038/nbt.3988>
    2. Traag, V. A., Waltman. L., van Eck, N. J., 2019.
       From Louvain to Leiden: guaranteeing well-connected communities.
       Scientific Reports 9, 5233.
       <https://doi.org/10.1038/s41598-019-41695-z>
    3. Csárdi, G., Nepusz, T., 2006.
       The igraph software package for complex network research.
       InterJournal Complex Systems, 1695.
    """
    check_executable_exists("mmseqs", "MMseqs2")

    if preset == PRESET.NT_PRECLUSTER:
        if min_seq_id is None:
            min_seq_id = 0.98
        if min_aln_len is None:
            min_aln_len = 100
        if seq_id_mode is None:
            seq_id_mode = 2
        if cov_mode is None:
            cov_mode = 0
        if coverage is None:
            coverage = 0.5
        if evalue is None:
            evalue = 1e-6
        if sensitivity is None:
            sensitivity = 7.5
        if cluster_mode is None:
            cluster_mode = 1
    elif preset == preset.NT_CLUSTER:
        check_executable_exists("seqkit")
        if sensitivity is None:
            sensitivity = 7.5
        if evalue is None:
            evalue = 1e-3
        if min_seq_id is None:
            min_seq_id = 0.40
        if max_seqs is None:
            max_seqs = 1_000_000
        if k is None:
            k = 5
    elif preset == PRESET.ORF:
        if min_seq_id is None:
            min_seq_id = 0.90
        if seq_id_mode is None:
            seq_id_mode = 2
        if cov_mode is None:
            cov_mode = 1
        if coverage is None:
            coverage = 0.333
        if evalue is None:
            evalue = 0.1
        if cluster_mode:
            cluster_mode = 2

    logging.info(f"Clustering{' with preset ' + preset if preset != 'none' else ''}...")

    # when the user hasn't specified a prefix, use the input FASTA file name
    if prefix is None:
        prefix = Path(fasta).stem

    base_command = f"easy-{'linclust' if lin else 'cluster'}"
    if preset == PRESET.NT_CLUSTER:
        base_command = f"easy-{'lin' if lin else ''}search"

    if preset == PRESET.NT_CLUSTER:
        # duplicate input sequences
        subprocess.run(
            f"seqkit concat --quiet {fasta} {fasta} > {Path(f'{fasta.stem}.duplicated.fasta')}",
            shell=True,
            check=True,
        )
        arg2 = str(fasta)
        arg3 = prefix + "_AvA"

    else:
        arg2 = prefix
        arg3 = ""

    logfile = Path(f"{prefix}.{preset + '.' if preset != 'none' else '' }mmseqs.txt")

    command = (
        f"mmseqs {base_command} "
        f"{'' if min_seq_id is None else f'--min-seq-id {min_seq_id} '}"
        f"{'' if min_aln_len is None else f'--min-aln-len {min_aln_len} '}"
        f"{'' if seq_id_mode is None else f'--seq-id-mode {seq_id_mode} '}"
        f"{'' if cov_mode is None else f'--cov-mode {cov_mode} '}"
        f"{'' if coverage is None else f'-c {coverage} '}"
        f"{'' if evalue is None else f'-e {evalue} '}"
        f"{'' if sensitivity is None or lin else f'-s {sensitivity} '}"
        f"{'' if cluster_mode is None else f'--cluster-mode {cluster_mode} '}"
        # nt-cluster is actually easy-search under the hood so it needs special params
        f"{'' if max_seqs is None else f'--max-seqs {max_seqs} '}"
        f"{'' if k is None else f'-k {k} '}"
        f"{'--search-type 3 ' if preset == PRESET.NT_CLUSTER else ''}"
        f"{'--format-output ' + NT_CLUSTER_COLNAMES + ' ' if preset == PRESET.NT_CLUSTER else ''}"
        # general options
        f"--threads {threads} "
        f"{fasta} {arg2} {arg3} {tmpdir} "
        f"> {logfile}"
    )
    try:
        logging.debug(f"{command=}")
        subprocess.run(
            command,
            shell=True,
            check=True,
        )

        # remove useless files
        Path(prefix + "_all_seqs.fasta").unlink(missing_ok=True)

        logging.done("Done clustering.")  # type: ignore
    except subprocess.CalledProcessError as e:
        logging.error(f"MMseqs failed with exit code {e.returncode}")
        if e.output:
            logging.error(f"MMseqs output: '{e.output.decode('utf-8').rstrip()}'")
        raise click.ClickException("MMseqs failed. See logs for error messages.")
    finally:
        # remove temporary files no matter what
        if tmpdir.exists():
            shutil.rmtree(tmpdir)

    # we might need to postprocess the output to make it consistent
    if preset == PRESET.NT_CLUSTER:
        AvA2cluster(path=Path(arg3), outfile=Path(prefix + "_cluster.tsv"))
        # we duplicated the input file so let's clean up
        # double check the name just to be sure
        if ".duplicated.fasta" in fasta.stem:
            fasta.unlink()
