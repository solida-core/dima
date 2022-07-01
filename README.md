# Snakemake workflow: DiMA
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.15.0-brightgreen.svg)](https://snakemake.bitbucket.io)

This workflow performs mapping of single-end and paired-end reads in fastq format against a reference genome to produce a deduplicated and recalibrated BAM file.

DiMA is part of the Snakemake-based pipelines collection [solida-core](https://github.com/solida-core) developed and manteined at [CRS4](https://www.crs4.it). 
<p align="center">
<img align="center" src="https://www.crs4.it/wp-content/uploads/2020/11/CRS4-1.jpg" width="200" height="80" alt="www.crs4.it"/>
</p>

## Authors

* Matteo Massidda (@massiddaMT)
* Rossano Atzeni (@ratzeni)

## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=solida-core/dima).

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

## INSTRUCTIONS
Create a virtual environment with the command:
```commandline
mamba create -c bioconda -c conda-forge --name snakemake snakemake=6.15 snakedeploy
```
and activate it:
```commandline
conda activate snakemake
```
We get some public data to test the pipeline. You can directly clone in this folder from github, just type:
```commandline
git clone https://github.com/solida-core/test-data-DNA.git
```
You can then perform the pipeline deploy defining a directory `my_dest_dir` for analysis output and a pipeline tag for a specific version:
```bash
snakedeploy deploy-workflow https://github.com/solida-core/dima 
                    my_desd_dir 
```
To run the pipeline, go inside the deployed pipeline folder and use the command:
```bash
snakemake --use-conda -p --cores all
```
You can generate analysis report with the command:
```bash
snakemake --report report.zip --cores all
```
