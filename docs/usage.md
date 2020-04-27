# hyperclust: Usage

## Input

### Sample input

The sample vcfs have to be provided through a index (csv) file which contains
the sample names, the vcf path, the path to the cnaFile information and the
purity of the sample.

Each column contains the path to:

- `sampleId` the name of the sample
- `vcfFile` path to the vcf file (only one sample, somatic)
- `cnaFile` a file with copy number alteration data
- `purity` a numeric value which corresponds to the tumor sample purity

The CNA file should be a `.tsv` file with the columns: chromosome, start, end,
total_cn, major_cn, minor_cn, star. The file can be obtained with tools like
[ASCAT](https://github.com/Crick-CancerGenomics/ascat).

As example, the index file should look like this.

```text
sampleId,vcfFile,cnaFile,purity
sID123,path/to/file.vcf.gz,path/to/file.cna.txt,0.75
```

### Fasta genome

The genome information needs to be provided through a fasta file. This file
can be hard masked (aka masked positions with N) to avoid randomization
of mutations into those regions.

The genome fasta file should be provided as a unique fasta file.

## Parameters

### Nextflow common params

Nextflow comes with a set of params to control some execution options, from
them, the more relevants are:

| Params |  description |
|----------|-------------------------|
| -w | Defines where the temporal files of the pipeline will be stored. |
| -c | Defines which configuration files should be used for the current execution |
| -profile | Defines which configuration profile should be used for the current execution |
| -resume | If to resume previously executed pipeline, intermediate files will be kept |
| -N | Email address to send execution report when the pipeline finishes or crashes, needs to be [configured in the system](https://www.nextflow.io/docs/latest/mail.html)|

### Pipeline specific params

The parameters which should be used in hyperclust are defined
below:

| Params | default | description |
|----------|---------|------------------------|
| --index | --- | Index file with vcf files and the corresponding sample names (csv format). |
| --genome| --- | Masked fasta file to use as reference sequence. |
| --dataset | "hartwig" | Different parsing options for common available datasets. (hartwig, pcawg_sanger and tcga_strelka) |
| --stratification | true | If stratify the samples according to clonality, strand and pair. (only true available rn) |
| --batchsize | int | [randommut](https://github.com/davidmasp/randommut) specific option. |
|  --intraBS |  int | [randommut](https://github.com/davidmasp/randommut) specific option. |
|  --ws  | int |      [randommut](https://github.com/davidmasp/randommut) specific option. |
|  --times | int |    [randommut](https://github.com/davidmasp/randommut) specific option. |
|  --clonal | bool |  [clustmut](https://github.com/davidmasp/clustMut) specific option. |

## Example

An example of a pipeline call would be:

```bash
nextflow run davidmasp/hyperclust \
    -profile slurm,conda \
    -w /path/to/scratch \
    -resume \
    -N you@gmail.com \
    --index samples.csv \
    --genome hg19.fa
```
