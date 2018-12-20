def get_fastq(wildcards,units,read_pair='fq1'):
    return units.loc[wildcards.unit,
                     [read_pair]].dropna()[0]


rule pre_rename_fastq_pe:
    input:
        r1=lambda wildcards: get_fastq(wildcards, units, 'fq1'),
        r2=lambda wildcards: get_fastq(wildcards, units, 'fq2')
    output:
        r1=temp(dima_path+"reads/{unit}-R1.fq.gz"),
        r2=temp(dima_path+"reads/{unit}-R2.fq.gz")
    shell:
        "ln -s {input.r1} {output.r1} &&"
        "ln -s {input.r2} {output.r2} "


rule trim_galore_pe:
    input:
        [dima_path+"reads/{unit}-R1.fq.gz", dima_path+"reads/{unit}-R2.fq.gz"]
    output:
        temp(dima_path+"reads/trimmed/{unit}-R1_val_1.fq.gz"),
        dima_path+"reads/trimmed/{unit}-R1.fq.gz_trimming_report.txt",
        temp(dima_path+"reads/trimmed/{unit}-R2_val_2.fq.gz"),
        dima_path+"reads/trimmed/{unit}-R2.fq.gz_trimming_report.txt"
    params:
        extra=config.get("rules").get("trim_galore_pe").get("arguments")
    log:
        dima_path+"logs/trim_galore/{unit}.log"
    benchmark:
        dima_path+"benchmarks/trim_galore/{unit}.txt"
    wrapper:
        "0.27.0/bio/trim_galore/pe"


rule post_rename_fastq_pe:
    input:
        r1=dima_path+"reads/trimmed/{unit}-R1_val_1.fq.gz",
        r2=dima_path+"reads/trimmed/{unit}-R2_val_2.fq.gz"
    output:
        r1=temp(dima_path+"reads/trimmed/{unit}-R1-trimmed.fq.gz"),
        r2=temp(dima_path+"reads/trimmed/{unit}-R2-trimmed.fq.gz")
    shell:
        "mv {input.r1} {output.r1} &&"
        "mv {input.r2} {output.r2} "

