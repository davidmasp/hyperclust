# hyperclust: Output

## Output directory

The output directory is the folder where all the output files will be stored.
It can be customized in the pipeline call as `--outdir outdir`.
Inside the output dir there will be also a tracedir which will store
execution information for the pipeline run and a DAG.

## Output files

The output directory of hyperclust will be structured in
different folders containing different output types.

| Folder | File Type | Description |
|---------|-----------|-----------------------|
| clustmut | `.rds` | Contains cluster calls from [clustmut](https://github.com/davidmasp/clustMut) |
| randommut | `.tsv` | Contains randomized positions per sample from [randommut](https://github.com/davidmasp/randommut) |
| stratification | various | Contains both plots from the clonal stratification and the stratification files |

## TMP folder

Some pipelines in development state will include a `TMP` folder which contains
intermediate files for debugging.
