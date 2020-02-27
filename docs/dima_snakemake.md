# DiMA
**DiMA** (DNA Mapping) is a pipeline for Next-Generation Sequencing data alignment.
_______________


## Requirements
DiMA pipeline contains the instruction to automatically download and install all required dependecies. 
The only required software is asked to install is [Conda](https://docs.conda.io/en/latest/miniconda.html).
 
## Usage
#### Preliminary Steps
For pipeline manual execution, just clone the git repository:
```bash
git clone https://github.com/solida-core/dima
```
The `dima/` folder contains all files required for pipeline execution, apart data-related files that the user have to include in the folder.

As mentioned above, Snakemake will manage all pipeline requirements, but we need to work in a virtual environment in which is installed Snakemake. The instructions are included in the `environment.yaml` file.

* ##### Create the project virtual environment:
```bash
conda env create --name my_environment_name --file environment.yaml
```
* ##### Activate the project virtual environment:
```bash
source activate my_environment_name
```
#### Run Snakemake

Now that the virtual environment is created and activated, it is possible to execute the pipeline.
Starting from final output indicated in the `rule_all` section of the [Snakefile](https://github.com/solida-core/docs/blob/master/pipeline_structure.md#snakefile), Snakemake will check for rules and input required for producing those files.
In this way, it is possible to perform a test-check called `dryrun`:
```bash
snakemake --snakefile Snakefile --configfile my_config.yaml --dryrun -d my_analysis_folder
```
The process will display, in absence of errors, all steps that need to be executed to obtain output indicated in the Snakefile's rule_all and the number of each rule execution repeat based on the number of samples.

The built graph can also be plotted with the following command:
```bash
snakemake --configfile my_config.yaml --dag | dot -Tpng > my_graph.png
```
Finally, the **execution of the pipeline** can be launched with the command:
```bash
snakemake --snakefile Snakefile --configfile my_config.yaml -d my_analysis_folder
```
For further [Snakemake options](https://snakemake.readthedocs.io/en/stable/api_reference/snakemake.html) type:
```bash
snakemake --help
```

