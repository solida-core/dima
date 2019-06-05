## Reference Files

[DiMA]() reference files are stored within a particular folder-tree common with all solida-core pipelines.

This organization must be kept because the pipelines have built-in functions which search for reference and accessory files in these folders.

This is a representation of it:
```
REFERENCE FOLDER
    │
    └── PROVIDER
            |
            └── RELEASE
                   |
                   ├── genome.fasta
                   |
                   ├── genome.idx
                   |    
                   └── known_varians
                            |
                            └── accessory files
```
There is a `main reference folder` in which can be present reference files from different `PROVIDER` (like ucsc or ncbi).

For each provider can be present multiple `RELEASES` of human genome, where are placed the relative files.

We refer to this first described three elemets as `REFERENCE`, `PROVIDER` and `RELEASE` and must be declared in the pipeline `config.yaml` file. 
Starting from this, a pipeline function is able to generally manage and use reference files from the correct provider and release.


The `RELEASE` folder is the directory in which are stored the reference genome files with relative indexes and, inside the "known_variants" folder, several GATK accessory files (useful for Variant Recalibration step).


## REFERENCE DOWNLOADER

To facilitate users in reference organization, we developed a simple script aimed at downloading and correctly placing reference files.
The script is part of a *[TOOLKIT](https://github.com/solida-core/toolkit#solida-core-toolkit)* where useful pipeline-related scripts will be collected.

**INSTALL and RUN:**  **[reference_organizer](https://github.com/solida-core/toolkit#reference_organizerpy)**

