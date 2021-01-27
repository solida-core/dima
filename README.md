[![depends](https://img.shields.io/badge/depends%20from-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![snakemake](https://img.shields.io/badge/snakemake-5.3-brightgreen.svg)](https://snakemake.readthedocs.io/en/stable/)
[![Travis](https://travis-ci.com/solida-core/dima.svg?branch=master)](https://travis-ci.com/solida-core/dima.svg?branch=master)

# DiMA
**DiMA** (DNA Mapping) is a pipeline for Next-Generation Sequencing data alignment.

All **[solida-core](https://github.com/solida-core)** workflows follow GATK Best Practices for Germline Variant Discovery, with the incorporation of further improvements and refinements after their testing with real data in various [CRS4 Next Generation Sequencing Core Facility](http://next.crs4.it) research sequencing projects.

Pipelines are based on [Snakemake](https://snakemake.readthedocs.io/en/stable/), a workflow management system that provides all the features needed to create reproducible and scalable data analyses.

Software dependencies are specified into the `environment.yaml` file and directly managed by Snakemake using [Conda](https://docs.conda.io/en/latest/miniconda.html), ensuring the reproducibility of the workflow on a great number of different computing environments such as workstations, clusters and cloud environments.


### Pipeline Overview
The pipeline workflow is aimed at [_Mapping_](docs/dima_workflow.md#mapping) paired-end reads in fastq format against a reference genome to produce a deduplicated and recalibrated BAM file.

Obtained BAM files can be then included in Variant Calling processes or visualized with tools like [IGV]().
DiMA pipeline is included in other solida-core pipelines that requires the mapping step (i.e. [DiVA]()). 
The standalone usage is recommended for analysis which requires BAM files and not a Variant Calling step.

A complete view of the analysis workflow is provided by the pipeline's [graph](images/dima.png).



### Pipeline Handbook
**DiMA** pipeline documentation can be found in the `docs/` directory:


1. [Pipeline Structure:](https://github.com/solida-core/docs/blob/master/pipeline_structure.md)
    * [Snakefile](https://github.com/solida-core/docs/blob/master/pipeline_structure.md#snakefile)
    * [Configfile](https://github.com/solida-core/docs/blob/master/pipeline_structure.md#configfile)
    * [Rules](https://github.com/solida-core/docs/blob/master/pipeline_structure.md#rules)
    * [Envs](https://github.com/solida-core/docs/blob/master/pipeline_structure.md#envs)
2. [Pipeline Workflow](docs/dima_workflow.md)
3. [Required Files:]()
    * [Reference files](docs/reference_files.md)
    * [User files](docs/user_files.md)
4. [Running the pipeline:]()
    * [Manual Snakemake Usage](docs/dima_snakemake.md)
    * [SOLIDA:]()
        * [CLI - Command Line Interface](https://github.com/solida-core/solida/blob/master/README.md)
        * [GUI - Graphical User Interface]()






### Contact us
[support@solida-core](mailto:m.massidda@crs4.it) 
