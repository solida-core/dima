rule samtools_sort:
   input:
       resolve_results_filepath(config.get("paths").get("results_dir"),"reads/aligned/{unit}_fixmate.cram")
   output:
       temp(resolve_results_filepath(config.get("paths").get("results_dir"),"reads/sorted/{unit}_sorted.cram"))
   conda:
       resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
   params:
       tmp_dir=tmp_path(path=config.get("paths").get("tmp_dir")),
       genome=config.get("resources").get("reference"),
       output_fmt="CRAM"
   benchmark:
       resolve_results_filepath(config.get("paths").get("results_dir"),"benchmarks/samtools/sort/{unit}.txt")
   threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
   shell:
       "samtools sort "
       "--threads {threads} "
       "-T {params.tmp_dir} "
       "-O {params.output_fmt} "
       "--reference {params.genome} "
       "-o {output} "
       "{input} "



rule samtools_merge:
    """
    Merge cram files for multiple units into one for the given sample.
    If the sample has only one unit, files will be copied.
    """
    input:
        lambda wildcards: get_units_by_sample(wildcards, samples,
                                              prefix=resolve_results_filepath(config.get("paths").get("results_dir"),'reads/sorted/'),
                                              suffix='_sorted.cram')
    output:
        temp(resolve_results_filepath(config.get("paths").get("results_dir"),"reads/merged/{sample}.cram"))
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    benchmark:
        resolve_results_filepath(config.get("paths").get("results_dir"),"benchmarks/samtools/merge/{sample}.txt")
    params:
        cmd='samtools',
        genome=config.get("resources").get("reference"),
        output_fmt="CRAM"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    script:
        "scripts/samtools_merge.py"


rule samtools_cram_to_bam:
    input:
        resolve_results_filepath(config.get("paths").get("results_dir"),"reads/merged/{sample}.cram")
    output:
        resolve_results_filepath(config.get("paths").get("results_dir"),"reads/merged/{sample}.bam")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    benchmark:
        resolve_results_filepath(config.get("paths").get("results_dir"),"benchmarks/samtools/cram_to_bam/{sample}.txt")
    params:
        genome=config.get("resources").get("reference"),
        output_fmt="BAM"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "samtools view -b "
        "--threads {threads} "
        "-T {params.genome} "
        "-o {output} "
        "-O {params.output_fmt} "
        "{input} "


rule samtools_index:
    input:
        resolve_results_filepath(config.get("paths").get("results_dir"),"reads/recalibrated/{sample}.dedup.recal.bam")
    output:
         resolve_results_filepath(config.get("paths").get("results_dir"),"reads/recalibrated/{sample}.dedup.recal.bam.bai")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/samtools.yaml")
    benchmark:
        resolve_results_filepath(config.get("paths").get("results_dir"),"benchmarks/samtools/index_2/{sample}.txt")
    shell:
        "samtools index "
        "{input} "
