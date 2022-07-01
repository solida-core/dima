
rule pre_rename_fastq_se:
    input:
        lambda wildcards: get_fastq(wildcards,units)
    output:
        r1=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/untrimmed/se/{unit}-R1.fq.gz")
    log:
        resolve_results_filepath(config.get("paths").get("results_dir"),"logs/bash/{unit}_cp.log")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/bash.yaml")
    shell:
        "cp {input} {output.r1} "
        ">& {log} "


rule trim_galore_se:
    input:
        rules.pre_rename_fastq_se.output
    output:
        fq1=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/trimmed/se/{unit}-R1_trimmed.fq.gz"),
        report_fq1=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/trimmed/se/{unit}-R1.fq.gz_trimming_report.txt")
    params:
        extra=config.get("params").get("trim_galore_se").get("arguments"),
        outdir=lambda w, output: os.path.dirname(output[0]),
        qc_dir=resolve_results_filepath(config.get("paths").get("results_dir"),"qc/fastqc")
        # outdir=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/trimmed/se/")
    log:
        resolve_results_filepath(config.get("paths").get("results_dir"),"logs/trim_galore/{unit}.log")
    benchmark:
        resolve_results_filepath(config.get("paths").get("results_dir"),"benchmarks/trim_galore/{unit}.txt")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/trim_galore.yaml")
    threads: (conservative_cpu_count(reserve_cores=2, max_cores=99))/2 if (conservative_cpu_count(reserve_cores=2, max_cores=99)) >2 else 1
    shell:
        "mkdir -p {params.qc_dir}; "
        "trim_galore "
        "{params.extra} "
        "--cores {threads} "
        "-o {params.outdir} "
        "{input} "
        ">& {log}"


rule post_rename_fastq_se:
    input:
        rules.trim_galore_se.output.fq1
    output:
        r1=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/se/trimmed/{unit}-R1-trimmed.fq.gz")
    log:
        resolve_results_filepath(config.get("paths").get("results_dir"),"logs/bash/{unit}_mv.log")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/bash.yaml")
    shell:
        "mv {input} {output.r1}"
        ">& {log} "