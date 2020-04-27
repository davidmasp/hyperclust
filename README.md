# hyperClust

**WIP WARNING**

## Overview

hyperClust is a pipeline for detection of somatic local hyper-mutation
(mutation clusters).
In brief, mutations that are located closer than expected to their
nearest neighbor (closest adjacent mutation) are considered clusters
of mutagenesis. These clusters can have different characteristics and
be generated from a sort of different mutational processes.

The pipeline is built in [nextflow](https://www.nextflow.io/),
a pipeline framework for reproducible research. The software handles
parallelization in multiple HPC scheduler or multicore machine.
See [installation](docs/install.md) for more details.

## Installation

### Requirements

A Unix system plus a nextflow installation
is required, other required software is listed here.

* [randommut](https://github.com/davidmasp/randommut)
* [clustMut](https://github.com/davidmasp/clustMut)

Installation instructions for these tools are available in each site.

## Usage

1. Install requirements
2. Build index files
3. Mask the reference genome
4. Run the pipeline

```bash
nextflow run davidmasp/hyperclust \
    -profile slurm,conda \
    -w /path/to/scratch \
    -resume \
    -N you@gmail.com \
    --index index.csv \
    --genome refseq.fa
```

## Documentation

A brief documentation can be found in the [`docs/`](docs) directory.

- [Installation](docs/install.md)
- [Configuration](docs/configuration.md)
- [Usage](docs/usage.md)
- [Output](docs/output.md)
- [Troubleshooting](docs/troubleshooting.md)

Please, if you find any issue or question report it in the
[issue section](https://github.com/davidmasp/hyperclust/issues).

## Credit

Author: David Mas-Ponte

