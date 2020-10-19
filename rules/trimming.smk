
def get_fastq(wildcards,units):
    # print(wildcards.unit)
    if units.loc[wildcards.unit,["fq2"]].isna().all():
        print("SE")
        # print(units.loc[wildcards.unit,["fq1"]].dropna()[0])
        return units.loc[wildcards.unit,["fq1"]].dropna()[0]
    else:
        print("PE")
        # print(units.loc[wildcards.unit,["fq1"]].dropna()[0],units.loc[wildcards.unit,["fq2"]].dropna()[0])
        return units.loc[wildcards.unit,["fq1"]].dropna()[0],units.loc[wildcards.unit,["fq2"]].dropna()[0]


rule pre_rename_fastq_pe:
    input:
        lambda wildcards: get_fastq(wildcards,units)
    output:
        r1="reads/untrimmed/{unit}-R1.fq.gz",
        r2="reads/untrimmed/{unit}-R2.fq.gz"
    shell:
        "ln -s {input[0]} {output.r1} &&"
        "ln -s {input[1]} {output.r2} "

rule pre_rename_fastq_se:
    input:
       lambda wildcards: get_fastq(wildcards,units)
    output:
        r1="reads/untrimmed/se/{unit}-R1.fq.gz"
    shell:
        "ln -s {input} {output.r1} "


rule trim_galore_pe:
    input:
        rules.pre_rename_fastq_pe.output
    output:
        temp("reads/trimmed/{unit}-R1_val_1.fq.gz"),
        "reads/trimmed/{unit}-R1.fq.gz_trimming_report.txt",
        temp("reads/trimmed/{unit}-R2_val_2.fq.gz"),
        "reads/trimmed/{unit}-R2.fq.gz_trimming_report.txt"
    params:
        extra=config.get("rules").get("trim_galore_pe").get("arguments"),
        outdir="reads/trimmed/"
    log:
        "logs/trim_galore/{unit}.log"
    benchmark:
        "benchmarks/trim_galore/{unit}.txt"
    conda:
        "../envs/trim_galore.yaml"
    threads: (conservative_cpu_count(reserve_cores=2, max_cores=99))/4
    shell:
        "mkdir -p qc/fastqc; "
        "trim_galore "
        "{params.extra} "
        "--cores {threads} "
        "-o {params.outdir} "
        "{input} "
        ">& {log}"


rule trim_galore_se:
    input:
        rules.pre_rename_fastq_se.output
    output:
        temp("reads/trimmed/se/{unit}-R1_trimmed.fq.gz"),
        "reads/trimmed/se/{unit}-R1.fq.gz_trimming_report.txt"
    params:
        extra=config.get("rules").get("trim_galore_se").get("arguments"),
        outdir="reads/trimmed/se/"
    log:
        "logs/trim_galore/{unit}.log"
    benchmark:
        "benchmarks/trim_galore/{unit}.txt"
    conda:
        "../envs/trim_galore.yaml"
    threads: (conservative_cpu_count(reserve_cores=2, max_cores=99))/2
    shell:
        "mkdir -p qc/fastqc; "
        "trim_galore "
        "{params.extra} "
        "--cores {threads} "
        "-o {params.outdir} "
        "{input} "
        ">& {log}"




rule post_rename_fastq_pe:
    input:
        rules.trim_galore_pe.output
    output:
        r1=temp("reads/trimmed/{unit}-R1-trimmed.fq.gz"),
        r2=temp("reads/trimmed/{unit}-R2-trimmed.fq.gz")
    shell:
        "mv {input[0]} {output.r1} &&"
        "mv {input[2]} {output.r2} "


rule post_rename_fastq_se:
    input:
        rules.trim_galore_se.output
    output:
        r1=temp("reads/se/trimmed/{unit}-R1-trimmed.fq.gz")
    shell:
        "mv {input[0]} {output.r1}"



def get_trimmed_reads(wildcards,units):
    print(wildcards.unit)
    if units.loc[wildcards.unit,["fq2"]].isna().all():
        # SE
        return rules.post_rename_fastq_se.output
    # PE
    else:
        return rules.post_rename_fastq_pe.output.r1,rules.post_rename_fastq_pe.output.r2
