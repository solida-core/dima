
rule pre_rename_fastq_pe:
    input:
        lambda wildcards: get_fastq(wildcards,units)
    output:
        r1=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/untrimmed/{unit}-R1.fq.gz"),
        r2=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/untrimmed/{unit}-R2.fq.gz")
    shell:
        "cp {input[0]} {output.r1} &&"
        "cp {input[1]} {output.r2} "


rule trim_galore_pe:
    input:
        rules.pre_rename_fastq_pe.output.r1,
        rules.pre_rename_fastq_pe.output.r2
    output:
        resolve_results_filepath(config.get("paths").get("results_dir"),"reads/trimmed/{unit}-R1_val_1.fq.gz"),
        resolve_results_filepath(config.get("paths").get("results_dir"),"reads/trimmed/{unit}-R1.fq.gz_trimming_report.txt"),
        resolve_results_filepath(config.get("paths").get("results_dir"),"reads/trimmed/{unit}-R2_val_2.fq.gz"),
        resolve_results_filepath(config.get("paths").get("results_dir"),"reads/trimmed/{unit}-R2.fq.gz_trimming_report.txt")
    params:
        extra=config.get("params").get("trim_galore_pe").get("arguments"),
        outdir=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/trimmed/")
    log:
        resolve_results_filepath(config.get("paths").get("results_dir"),"logs/trim_galore/{unit}.log")
    benchmark:
        resolve_results_filepath(config.get("paths").get("results_dir"),"benchmarks/trim_galore/{unit}.txt")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/trim_galore.yaml")
    threads: (conservative_cpu_count(reserve_cores=2, max_cores=99))/4 if (conservative_cpu_count(reserve_cores=2, max_cores=99)) >4 else 1
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
        r1=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/trimmed/{unit}-R1-trimmed.fq.gz"),
        r2=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/trimmed/{unit}-R2-trimmed.fq.gz")
    shell:
        "mv {input[0]} {output.r1} &&"
        "mv {input[2]} {output.r2} "