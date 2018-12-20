

rule mark_duplicates:
    input:
        dima_path+"reads/merged/{sample}.bam"
    output:
        bam=temp(dima_path+"reads/dedup/{sample}.dedup.bam"),
        metrics=dima_path+"reads/dedup/{sample}.metrics.txt"
    log:
        dima_path+"logs/picard/MarkDuplicates/{sample}.log"
    benchmark:
        dima_path+"benchmarks/picard/MarkDuplicates/{sample}.txt"
    params:
        java_params(tmp_dir=config.get("paths").get("to_tmp"), multiply_by=5),
        config.get("rules").get("picard_MarkDuplicates").get("arguments"),
        lambda wildcards: get_odp(wildcards, samples, 'odp')
    wrapper:
        "0.27.0/bio/picard/markduplicates"