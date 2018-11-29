import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.1.2")

##### load config and sample sheets #####

#configfile: "config.yaml"

samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col=["unit"], dtype=str)
#units.index = units.index.set_levels([i.astype(str) for i in units.index.levels]) # enforce str in index


##### local rules #####

localrules: all, pre_rename_fastq_pe, post_rename_fastq_pe


##### target rules #####

rule all:
    input:
        expand("reads/merged/{sample.sample}.cram.crai",
              sample=samples.reset_index().itertuples()),
#        expand("reads/trimmed/{unit.unit}-R{read}-trimmed.fq.gz",
#               unit=units.reset_index().itertuples(),
#               read=[1, 2]),
#        expand("reads/aligned/{unit.unit}_fixmate.cram",
#               unit=units.reset_index().itertuples()),
#        expand("reads/sorted/{unit.unit}_sorted.cram",
#               unit=units.reset_index().itertuples()),
#        expand("reads/merged/{sample.sample}.cram",
#               sample=samples.reset_index().itertuples()),
#        expand("reads/dedup/{sample.sample}.dedup.bam",
#            sample=samples.reset_index().itertuples()),
        expand("reads/recalibrated/{sample.sample}.dedup.recal.bam",
            sample=samples.reset_index().itertuples()),



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

