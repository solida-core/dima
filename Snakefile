import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")

##### load config and sample sheets #####

#configfile: "config.yaml"

samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col=["unit"], dtype=str)

##### local rules #####

localrules: all, pre_rename_fastq_pe, post_rename_fastq_pe


##### target rules #####

rule all:
    input:
        expand("reads/recalibrated/{sample.sample}.dedup.recal.bam",
            sample=samples.reset_index().itertuples()),
        "check_files.done"



##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3:4.4.10"


##### load rules #####

include_prefix="rules"

include:
    include_prefix + "/functions.py"
include:
    include_prefix + "/trimming.smk"
include:
    include_prefix + "/alignment.smk"
include:
    include_prefix + "/samtools.smk"
include:
    include_prefix + "/picard.smk"
include:
    include_prefix + "/bsqr.smk"


dima_path = ""

