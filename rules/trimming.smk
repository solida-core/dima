def get_fastq(wildcards,units,read_pair='fq1'):
    return units.loc[wildcards.unit,
                     [read_pair]].dropna()[0]


rule pre_rename_fastq_pe:
    input:
        r1=lambda wildcards: get_fastq(wildcards, units_pe, 'fq1'),
        r2=lambda wildcards: get_fastq(wildcards, units_pe, 'fq2')
    output:
        r1="reads/untrimmed/{unit}-R1.fq.gz",
        r2="reads/untrimmed/{unit}-R2.fq.gz"
    shell:
        "ln -s {input.r1} {output.r1} &&"
        "ln -s {input.r2} {output.r2} "

rule pre_rename_fastq_se:
    input:
        r1=lambda wildcards: get_fastq(wildcards, units_se, 'fq1')
    output:
        r1="reads/untrimmed/{unit}-R1.fq.gz"
    shell:
        "ln -s {input.r1} {output.r1} "


rule trim_galore_pe:
    input:
        ["reads/untrimmed/{unit}-R1.fq.gz", "reads/untrimmed/{unit}-R2.fq.gz"]
    output:
        temp("reads/trimmed/{unit}-R1_val_1.fq.gz"),
        "reads/trimmed/{unit}-R1.fq.gz_trimming_report.txt",
        temp("reads/trimmed/{unit}-R2_val_2.fq.gz"),
        "reads/trimmed/{unit}-R2.fq.gz_trimming_report.txt"
    params:
        extra=config.get("rules").get("trim_galore_pe").get("arguments")
    log:
        "logs/trim_galore/{unit}.log"
    benchmark:
        "benchmarks/trim_galore/{unit}.txt"
    wrapper:
        config.get("wrappers").get("trim_galore_pe")

rule trim_galore_se:
    input:
        "reads/untrimmed/merged/{unit}-R1.fq.gz"
    output:
        temp("reads/trimmed/{unit}_trimmed.fq.gz"),
        "reads/trimmed/{unit}.fq.gz_trimming_report.txt"
    params:
        extra=config.get("rules").get("trim_galore_se").get("arguments")
    log:
        "logs/trim_galore/{unit}.log"
    benchmark:
        "benchmarks/trim_galore/{units}.txt"
    wrapper:
        config.get("wrappers").get("trim_galore_se")



rule post_rename_fastq_pe:
    input:
        r1="reads/trimmed/{unit}-R1_val_1.fq.gz",
        r2="reads/trimmed/{unit}-R2_val_2.fq.gz"
    output:
        r1="reads/trimmed/{unit}-R1-trimmed.fq.gz",
        r2="reads/trimmed/{unit}-R2-trimmed.fq.gz"
    shell:
        "mv {input.r1} {output.r1} &&"
        "mv {input.r2} {output.r2} "

rule post_rename_fastq_se:
    input:
        r1="reads/trimmed/{sample}_trimmed.fq.gz"
    output:
        r1="reads/trimmed/{sample}-R1-trimmed.fq.gz"
    shell:
        "mv {input.r1} {output.r1}"