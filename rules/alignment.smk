
rule bwa_mem:
    input:
        lambda wildcards: get_trimmed_reads(wildcards,units)
    output:
        temp("reads/aligned/{unit}_fixmate.cram")
    conda:
        "../envs/bwa_mem.yaml"
    params:
        sample=lambda wildcards: '.'.join(wildcards.unit.split('.')[2:]),
        custom=config.get("rules").get("bwa-mem").get("arguments"),
        platform=config.get("rules").get("bwa-mem").get("platform"),
        platform_unit=lambda wildcards: '.'.join(wildcards.unit.split('.')[:-1]),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        output_fmt="CRAM"
    log:
        "logs/bwa_mem/{unit}.log"
    benchmark:
        "benchmarks/bwa/mem/{unit}.txt"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        'bwa mem {params.custom} '
        r'-R "@RG\tID:{wildcards.unit}\tSM:{params.sample}\tPL:{params.platform}\tLB:lib1\tPU:{params.platform_unit}" '
        '-t {threads} {params.genome} {input} 2> {log} '
        '|samtools fixmate --threads {threads} '
        '-O {params.output_fmt} '
        '--reference {params.genome} '
        '- {output} '
