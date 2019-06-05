# DiMA WORKFLOW
**DiMA** (DNA Mapping) is a pipeline for Next-Generation Sequencing data alignment.
______________________________
The pipeline workflow is aimed at [_Mapping_](docs/dima_workflow.md#mapping) paired-end reads in fastq format against a reference genome to produce a deduplicated and recalibrated BAM file.

Obtained BAM files can be then included in Variant Calling processes or visualized with tools like [IGV]().
DiMA pipeline is included in other solida-core pipelines that requires the mapping step (i.e. [DiVA]()). 
The standalone usage is recommended for analysis which requires BAM files and not a Variant Calling step.

A complete view of the analysis workflow is provided by the pipeline's [graph](images/dima.png).

_________________________________

### Mapping
* Trimming
* Alingment
* MarkDuplicates
* Base Recalibration


