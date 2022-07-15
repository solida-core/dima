
rule jvarkit_target_coverage:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/recalibrated/{sample}.dedup.recal.bam",
        ),
    output:
        tsv=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "results/stats/{sample}_target_coverage.tsv",
        ),
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=1),
        genome=config.get("resources").get("reference"),
        intervals=config.get("resources").get("bed"),
        mincov=config.get("params").get("jvarkit").get("min_coverage"),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    conda:
        "../envs/jvarkit.yaml"
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/jvarkit/{sample}_coverage.log",
        ),
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "bamstats04 "
        "-B {params.intervals} "
        "--cov {params.mincov} "
        "-R {params.genome} "
        "{input.bam} "
        "-o {output.tsv} "
        ">& {log}"


rule parse_jvarkit_coverage:
    input:
        target=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "results/stats/{sample}_target_coverage.tsv",
        ),
    output:
        target=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "results/tsv/{sample}.target_coverage.tsv",
        ),
        xlsx=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "results/xlsx/{sample}.coverage.xlsx",
        ),
    params:
        target_intervals=config.get("resources").get("bed"),
    conda:
        "../envs/r_env.yaml"
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/jvarkit/{sample}_parsing.log",
        ),
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    script:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/scripts/parse_results.R"
        )
