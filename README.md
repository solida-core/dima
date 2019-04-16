# DiMA 
Snakemake pipeline to map DNA datasets to a given reference genome using [BWA 
MEM](https://github.com/lh3/bwa/) and [Samtools](http://www.htslib.org/).  
Input datsets are fastq files  gzipped. Can be organized in 
multiple 
units for each sample or also only one unit per sample; in the last case, file will be 
copied.  
Output can be produced in BAM or CRAM (_default_) format.

## Workflow
![Dima dag](images/dima.png)

## Requirements
The pipeline's requirements are specified into the _environment.yml_ file and 
packages dependency are resolved using [Conda](https://conda.io/miniconda.html). 

## Usage

### Manual deployment (a.k.a. hard way)

Clone the repository and cd in it
```bash
git clone https://bitbucket.org/biopipelines/dima
cd dima
```

Edit the configuration file and the Snakefile to match your environment  
```
nano config.test.data.json   
nano Snakefile
```

Create conda environment  
``` 
conda env create -n dima --file environment.yaml
```

then activate it  
```
source activate dima
```

Launch Snakemake  
```
snakemake --use-conda --configfile config.test.data.json
```

### Automatic deployment (a.k.a. easy way)

Use [Solida](https://bitbucket.org/biopipelines/solida).

### Output

Default output of the pipeline in CRAM format, but can changed easily in BAM 
editing the Snakefile (look for OUTPUT_FORMAT variable) 


## Contributing

Contributions from everyone and anyone are welcome.  
Fork this repository, make your changes and create a Pull Request. 
Then one of the maintainers will review your changes.  
When all comments have been addressed and all tests pass, your changes will 
be merged.

