
rule pre_rename_fastq_pe:
    input:
        # lambda wildcards: get_fastq(wildcards,units)
        r1=lambda wildcards: get_a_fastq(wildcards, units, fq="fq1"),
        r2=lambda wildcards: get_a_fastq(wildcards, units, fq="fq2")
    output:
        r1=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/untrimmed/{unit}-R1.fq.gz"),
        r2=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/untrimmed/{unit}-R2.fq.gz")
    log:
        resolve_results_filepath(config.get("paths").get("results_dir"),"logs/bash/pe/{unit}_cp.log")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/bash.yaml")
    shell:
        "cp {input.r1} {output.r1} &&"
        "cp {input.r2} {output.r2} "
        ">& {log} "


rule trim_galore_pe:
    input:
        rules.pre_rename_fastq_pe.output.r1,
        rules.pre_rename_fastq_pe.output.r2
    output:
        fq1=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/trimmed/{unit}-R1_val_1.fq.gz"),
        report_r1=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/trimmed/{unit}-R1.fq.gz_trimming_report.txt"),
        fq2=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/trimmed/{unit}-R2_val_2.fq.gz"),
        report_r2=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/trimmed/{unit}-R2.fq.gz_trimming_report.txt")
    params:
        extra=config.get("params").get("trim_galore_pe").get("arguments"),
        outdir= lambda w,output: os.path.dirname(output[0]),
        qc_dir=resolve_results_filepath(config.get("paths").get("results_dir"),"qc/fastqc")
    log:
        resolve_results_filepath(config.get("paths").get("results_dir"),"logs/trim_galore/{unit}.log")
    benchmark:
        resolve_results_filepath(config.get("paths").get("results_dir"),"benchmarks/trim_galore/{unit}.txt")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/trim_galore.yaml")
    threads: (conservative_cpu_count(reserve_cores=2, max_cores=99))/4 if (conservative_cpu_count(reserve_cores=2, max_cores=99)) >4 else 1
    shell:
        "mkdir -p {params.qc_dir}; "
        "trim_galore "
        "{params.extra} "
        "--cores {threads} "
        "-o {params.outdir} "
        "{input} "
        ">& {log}"


rule post_rename_fastq_pe:
    input:
        r1=rules.trim_galore_pe.output.fq1,
        r2=rules.trim_galore_pe.output.fq2
    output:
        r1=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/trimmed/{unit}-R1-trimmed.fq.gz"),
        r2=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/trimmed/{unit}-R2-trimmed.fq.gz")
    log:
        resolve_results_filepath(config.get("paths").get("results_dir"),"logs/bash/pe/{unit}_mv.log")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/bash.yaml")
    shell:
        "mv {input.r1} {output.r1} &&"
        "mv {input.r2} {output.r2} "
        ">& {log} "