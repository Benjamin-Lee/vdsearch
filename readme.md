# `vdsearch`

Utilities for finding and analyzing viroid-like circular RNAs.

**Usage**:

```console
$ vdsearch [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--version`: Show the version and exit.
* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `dedup`: Remove duplicate sequences.
* `download-cms`: Download the latest covariance matrices for...
* `download-viroiddb`: Download the latest ViroidDB dataset.
* `easy-search`: Search viroid-like sequences.
* `find-circs`: Find circular RNA sequences.
* `infernal`: Run Infernal on the sequences.
* `mmseqs`: Run MMseqs on the sequences.
* `remove-rz-only-hits`: Remove sequences that only hit ribozymes.

## `vdsearch dedup`

Remove duplicate sequences.

**Usage**:

```console
$ vdsearch dedup [OPTIONS] FASTA
```

**Arguments**:

* `FASTA`: Path to .fasta or .fasta.gz file  [required]

**Options**:

* `--help`: Show this message and exit.

## `vdsearch download-cms`

Download the latest covariance matrices for ribozymes.

**Usage**:

```console
$ vdsearch download-cms [OPTIONS]
```

**Options**:

* `--help`: Show this message and exit.

## `vdsearch download-viroiddb`

Download the latest ViroidDB dataset.

**Usage**:

```console
$ vdsearch download-viroiddb [OPTIONS]
```

**Options**:

* `--help`: Show this message and exit.

## `vdsearch easy-search`

Search viroid-like sequences.

**Usage**:

```console
$ vdsearch easy-search [OPTIONS] FASTA
```

**Arguments**:

* `FASTA`: Path to .fasta or .fasta.gz file  [required]

**Options**:

* `--reference-db FILE`: Path to FASTA-formatted reference viroid database. If none is provided, the latest ViroidDB will be used.
* `--reference-cms FILE`: Path to Infernal ribozyme families. If none is provided, the latest ViroidDB will be used.
* `--help`: Show this message and exit.

## `vdsearch find-circs`

Find circular RNA sequences.

**Usage**:

```console
$ vdsearch find-circs [OPTIONS] FASTA OUTPUT
```

**Arguments**:

* `FASTA`: Path to .fasta or .fasta.gz file  [required]
* `OUTPUT`: Path to output file  [required]

**Options**:

* `--help`: Show this message and exit.

## `vdsearch infernal`

Run Infernal on the sequences.

**Usage**:

```console
$ vdsearch infernal [OPTIONS] FASTA
```

**Arguments**:

* `FASTA`: Path to .fasta or .fasta.gz file  [required]

**Options**:

* `--help`: Show this message and exit.

## `vdsearch mmseqs`

Run MMseqs on the sequences.

**Usage**:

```console
$ vdsearch mmseqs [OPTIONS] FASTA
```

**Arguments**:

* `FASTA`: Path to .fasta or .fasta.gz file  [required]

**Options**:

* `--help`: Show this message and exit.

## `vdsearch remove-rz-only-hits`

Remove sequences that only hit ribozymes.

**Usage**:

```console
$ vdsearch remove-rz-only-hits [OPTIONS] FASTA
```

**Arguments**:

* `FASTA`: Path to .fasta or .fasta.gz file  [required]

**Options**:

* `--help`: Show this message and exit.

