# vdsearch: find viroids and viroid-like circRNAs

Finding viroids used to be a difficult task. Now it can be done in one line of code:

```sh
vdsearch easy-search my-transcriptome.fasta
```

## Features

- End-to-end viroid-like RNA detection with one command
- Finds viroids (both families), satellite RNAs, ribozyviruses, and more
- _De novo_ and reference-free
- Works for both transcriptomes and metatranscriptomes

## Screenshots

![](screenshot.png)

## Example usage

- Find viroid-like cccRNAs in a transcriptome or metatranscriptome

  ```sh
  vdsearch easy-search my-transcriptome.fasta
  ```

- Detect circularity from assembled contigs

  ```sh
  vdsearch find-circs my-transcriptome.fasta
  ```

- Cluster circRNAs
  ```sh
  vdsearch cluster --preset nt-cluster circRNAs.fasta
  ```

## Installation

### Using Docker

The easiest way to use `vdsearch` is to use the [Docker image](https://hub.docker.com/r/benjaminlee/vdsearch). This will install all dependencies for you and run the code in a virtual environment.

```bash
docker pull benjamindlee/vdsearch
```

Note: When using Docker, you will have to mount the directory containing your input files to the container. For example, if you have a file called `my-transcriptome.fasta` in the current directory, you can run the following command:

```bash
docker run -v $(pwd):/data benjamindlee/vdsearch easy-search /data/my-transcriptome.fasta
```

You may have to configure Docker to use more memory and CPUs. For example, if you have 16 CPUs and 64 GB of RAM, you can run the following command:

```bash
docker run -v $(pwd):/data --cpus 16 -m 64G benjamindlee/vdsearch easy-search /data/my-transcriptome.fasta
```

### From source

If you want to run the code directly from source and will install dependencies (`seqkit`, etc.) yourself, please clone the repository and run the following commands:

```bash
git clone https://github.com/Benjamin-Lee/vdsearch.git
cd vdsearch

# activate a virtual environment (optional)
python3 -m venv venv
source venv/bin/activate

# install dependencies
pip install nimporter
pip install -e .

# test the installation
vdsearch --help
```

## Contributing

Contributions are always welcome!

See `contributing.md` for ways to get started.

## Support

For support, email [benjamin.lee@chch.ox.ac.uk](mailto:benjamin.lee@chch.ox.ac.uk) or [join our Slack channel](https://viroids.org/community).

## License

[MIT](https://choosealicense.com/licenses/mit/)
