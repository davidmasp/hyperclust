# HyperMut

HyperMut is a software for detection of somatic local hyper-mutation.
In brief, mutations that are located closer than expected to their
nearest neighbor (closest adjacent mutation) are considered clusters
of mutagenesis. These clusters can have different characteristics and
be generated from a sort of different mutational processes.

As overview of this you can check the recently published paper from our lab.

> Supek, F., & Lehner, B. (2017). Clustered mutation signatures reveal that error-prone DNA repair targets mutations to active genes. Cell, 170(3), 534-547.

## Installation

### Requirements

HyperMut is built with [nextflow](https://www.nextflow.io/) so an Unix system
plus a nextflow installation
is required, other required software can be listed here, however, installation
of the generated conda environtment by nextflow should be sufficient.

* randommut (python3 vXXX)
* clustMut (R v3.XXX)
* (optional[^1]) pyclone (python2)

[^1]: It is possible to run a self-contained clonal fraction estimator developed
in house

## Usage

## Documentation

## Output
