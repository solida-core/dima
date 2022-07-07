
rule bwa_mem:
    input:
        lambda wildcards: get_trimmed_reads(wildcards, units),
    output:
        resolve_results_filepath(
            config.get("paths").get("results_dir"), "reads/aligned/{unit}_fixmate.cram"
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/bwa_mem.yaml"
        )
    params:
        sample=lambda wildcards: ".".join(wildcards.unit.split(".")[2:]),
        custom=config.get("params").get("bwa-mem").get("arguments"),
        platform=config.get("params").get("bwa-mem").get("platform"),
        platform_unit=lambda wildcards: ".".join(wildcards.unit.split(".")[:-1]),
        genome=config.get("resources").get("reference"),
        output_fmt="CRAM",
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            resolve_results_filepath(
                config.get("paths").get("results_dir"), "logs/bwa_mem/{unit}.log"
            ),
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("results_dir"), "benchmarks/bwa/mem/{unit}.txt"
        )
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "bwa mem {params.custom} "
        r'-R "@RG\tID:{wildcards.unit}\tSM:{params.sample}\tPL:{params.platform}\tLB:lib1\tPU:{params.platform_unit}" '
        "-t {threads} {params.genome} {input} 2> {log} "
        "|samtools fixmate --threads {threads} "
        "-O {params.output_fmt} "
        "--reference {params.genome} "
        "- {output} "
