# hyperclust: Installation

## Nextflow

[Nextflow](https://www.nextflow.io/) is a pipeline software that permits
the execution of highly customizable and complex bioinformatic pipelines in
any HPC system and the cloud.

The hyperclust pipeline uses nextflow as pipeline framework so in
order to run this pipeline you will need to install the nextflow first.
To obtain updated guidelines of how to install nextflow follow visit
 [here](https://www.nextflow.io/).

Basic instructions:

```bash
# make sure java version is 1.8 or higher
java -version
cd PATH/in/your/bin
# this creates a nextflow executable
curl -s https://get.nextflow.io | bash
./nextflow run hello
```

The pipeline can be installed locally (the run directory will be the same
as the pipeline directory) or user-wide (the pipeline will be installed
in the home directory and it can be used system-wide)

## Local Installation

In order to use a local insallation you should clone this repository
and update the configuration options.

```bash
git clone https://github.com/davidmasp/hyperclust
```

This install option is suited if you need to adapt the configuration files
to your system.

If you choose to install the pipeline like this you should run it like
`nextflow run main.nf [OPT]`

### Use a given release

To obtain a given release use the tagged versions.

```bash
git checkout tags/<version> -b master
```

## Remote Run

You can also make a user-wise installation of the pipeline, or run it
remotely.

If pulling (install), nextflow will create a directory in your `$HOME` and
clone the pipeline.
If running (running), nextflow will do the same but check the installed version
is the latest. For more detailed info, see
[here](https://www.nextflow.io/docs/latest/sharing.html).

```bash
# for user-wide install
nextflow pull davidmasp/hyperclust

# for directly running the pipeline from github
nextflow run davidmasp/hyperclust [OPT]
```

If you choose to install the pipeline like this you should run it like
`nextflow run davidmasp/hyperclust [OPT]`
