

rule mark_duplicates:
    input:
        "reads/merged/{sample}.bam"
    output:
        bam=temp("reads/dedup/{sample}.dedup.bam"),
        bai=temp("reads/dedup/{sample}.dedup.bai"),
        metrics="reads/dedup/{sample}.metrics.txt"
    log:
        "logs/picard/MarkDuplicates/{sample}.log"
    conda:
        "../envs/picard.yaml"
    benchmark:
        "benchmarks/picard/MarkDuplicates/{sample}.txt"
    params:
        java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        config.get("rules").get("picard_MarkDuplicates").get("arguments"),
        lambda wildcards: get_odp(wildcards, samples, 'odp')
    shell:
        "picard MarkDuplicates "
        "{params} "
        "INPUT={input} "
        "OUTPUT={output.bam} "
        "METRICS_FILE={output.metrics} "
        ">& {log}"
