# Introduction
This pipeline computes synteny between two genomes. Synteny analysis pipeline implemented in a scalable, parallelized manner.
It uses mummer4 to compute synteny, but user could easily replace it with their own program.

## Input
This pipeline requires two genome files, query/focal genome and target genome as input.
Additionally, users can supply a list of query and target chromosomes to use. The synteny analysis will then be limited to these chromosomes.

## Output
The  pipeline is configured to output `.syn` files which could be directly used with the [fagin](https://github.com/lijing28101/fagin) program.


# Setup

Create the conda environment present in the `environment.yaml` file.

```
conda env create -f environment.yaml
```

## Configuration

The `config.yaml` files contains all the parameters for the pipeline. Please edit the `config.yaml` as required.

# Executing pipeline

## Using single node
To execute the pipeline on a single node use the following command:

```
snakemake -j 2 --use-conda --conda-frontend conda
```

Above command will execute the pipeline using 2 cores i.e. the transcript assembly step will run 2 samples in parallel.

## Scaling pipeline on HPC

To execute the pipeline on an HPC with mutiple nodes, execute as

```
snakemake -j 20 --profile snakemake_config/slurm
```

Above command will execute the pipeline ans schedule 20 jobs in parallel. 
The snakemake profile, `snakemake_config/slurm` could be modified to fine tune resource usage.
Above command is compatible for `SLURM` job-scheduler, you may need to change the parameters based on the system.
Specifically, edit the `snakemake_config/slurm/config.yaml` and `snakemake_config/slurm/cluster_config.yml`

