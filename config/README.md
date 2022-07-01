
# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample Files

Add samples to `config/samples.tsv` and `config/units.tsv`. 

## units
In the `units` file each row contain a single **unit** for a given sample, with fastq files for read1 and, if present, read2.
The file have 4 tab-separated columns: 
* `sample`: generic sample name, it will be reported also in the `samples.tsv` file. Can be present multiple times in different rows if the sample has multiple units.
* `unit`: is a unique identifier for the unit. A unit, i.e. `HSQ1008_141.L001.NA20804` is composed by 3 different parts:
  * `flowcell_id`: an id for the flowcell or instrument
  * `lane`: the lane id in which that unit was sequenced
  * `sample_id`: the sample id (the same in the first column of the file)
* `fq1`: absolute path of fastq file containing read1
* `fq2`: absolute path of fastq file containing read2, leave blank if SE.

An example ``config/units.tsv`` is reported below:
```
sample	unit	fq1	fq2
NA20804	HSQ1008_141.L001.NA20804	/abs_path/NA20804.L001.R1.fq.gz
NA20804	HSQ1008_141.L005.NA20804	/abs_path/NA20804.L005.R1.fq.gz	/abs_path/NA20804.L005.R2.fq.gz
```

## samples
In the `samples` file each row contain information for a single **sample**, with indication of all its units.
The file have 3 tab-separated columns: 
* `sample`: generic sample name, the same indicated in the `units.tsv` file. 
* `odp`: Optical Duplicate Distance used for Picard MarkDuplicates, depending on flowcell type.
* `units`: comma separated list of units (the unit name reported in the units file) for a given sample.

An example ``config/samples.tsv`` is reported below:
```
sample	odp	units
NA20804	100	HSQ1008_141.L001.NA20804,HSQ1008_141.L005.NA20804,HSQ1009_88.L007.NA20804
NA20806	100	HSQ1009_86.L001.NA20806
```


