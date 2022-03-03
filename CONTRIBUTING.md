## Adding new wrappers

If you add new wrappers, please try to mirror the original tool's command line flags.
However, there are a few exceptions:

1. Flags with underscores are not supported by Typer and are converted.
   For example, `--output_dir` is not supported even if the function takes `output_dir` as an argument.
2. The number of cores to be used must be specified with the `--threads` flag for consistency.
   Use the `Threads` object from `vdsearch.types` as the argument.

Be sure to include a citation for the original tool in the `References` section of the docstring.
The reference should be in Elsevier Harvard style.
