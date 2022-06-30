

rule mark_duplicates:
    input:
        resolve_results_filepath(config.get("paths").get("results_dir"),"reads/merged/{sample}.bam")
    output:
        bam=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/dedup/{sample}.dedup.bam"),
        bai=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/dedup/{sample}.dedup.bai"),
        metrics=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/dedup/{sample}.metrics.txt")
    log:
        resolve_results_filepath(config.get("paths").get("results_dir"),"logs/picard/MarkDuplicates/{sample}.log")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/picard.yaml")
    benchmark:
        resolve_results_filepath(config.get("paths").get("results_dir"),"benchmarks/picard/MarkDuplicates/{sample}.txt")
    params:
        java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=5),
        config.get("params").get("picard_MarkDuplicates").get("arguments"),
        lambda wildcards: get_odp(wildcards, samples, 'odp')
    shell:
        "picard MarkDuplicates "
        "{params} "
        "INPUT={input} "
        "OUTPUT={output.bam} "
        "METRICS_FILE={output.metrics} "
        ">& {log}"
