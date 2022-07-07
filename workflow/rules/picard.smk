

rule mark_duplicates:
    input:
        resolve_results_filepath(
            config.get("paths").get("results_dir"), "reads/merged/{sample}.bam"
        ),
    output:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"), "reads/dedup/{sample}.dedup.bam"
        ),
        bai=resolve_results_filepath(
            config.get("paths").get("results_dir"), "reads/dedup/{sample}.dedup.bai"
        ),
        metrics=resolve_results_filepath(
            config.get("paths").get("results_dir"), "reads/dedup/{sample}.metrics.txt"
        ),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/picard/MarkDuplicates/{sample}.log",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/picard.yaml"
        )
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "benchmarks/picard/MarkDuplicates/{sample}.txt",
        )
    params:
        java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=1),
        config.get("params").get("picard_MarkDuplicates").get("arguments"),
        lambda wildcards: get_odp(wildcards, samples, "odp"),
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "picard MarkDuplicates "
        "{params} "
        "INPUT={input} "
        "OUTPUT={output.bam} "
        "METRICS_FILE={output.metrics} "
        ">& {log}"


rule picard_pre_HsMetrics_probes:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/recalibrated/{sample}.dedup.recal.bam",
        ),
    output:
        probes=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "references/{sample}_probes_header",
            )
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/samtools.yaml"
        )
    params:
        probes=config.get("resources").get("probes"),
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/picard/preHsMetrics/{sample}.prehsmetrics_probes.log",
        ),
    shell:
        "samtools view -H  {input.bam} | cat - {params.probes} > {output.probes} "
        ">& {log} "


rule picard_pre_HsMetrics_target:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/recalibrated/{sample}.dedup.recal.bam",
        ),
    output:
        hsTarget=temp(
            resolve_results_filepath(
                config.get("paths").get("results_dir"),
                "references/{sample}_hsTarget_header",
            )
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/samtools.yaml"
        )
    params:
        hsTarget=config.get("resources").get("bed"),
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/picard/preHsMetrics/{sample}.prehsmetrics_target.log",
        ),
    shell:
        "samtools view -H  {input.bam} | cat - {params.hsTarget} > {output.hsTarget} "
        ">& {log} "


rule picard_HsMetrics:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/recalibrated/{sample}.dedup.recal.bam",
        ),
        probes=resolve_results_filepath(
            config.get("paths").get("results_dir"), "references/{sample}_probes_header"
        ),
        hsTarget=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "references/{sample}_hsTarget_header",
        ),
    output:
        metrics=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/recalibrated/{sample}.dedup.recal.hs.txt",
        ),
        target=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/recalibrated/{sample}.per_target_coverage.txt",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/picard.yaml"
        )
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=1),
        params=config.get("params").get("picard_HSmetrics").get("arguments"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/picard/HsMetrics/{sample}.hsmetrics.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "benchmarks/picard/HsMetrics/{sample}.txt",
        )
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "picard {params.custom} CollectHsMetrics "
        "{params.params} "
        "INPUT={input.bam} "
        "OUTPUT={output.metrics} "
        "BAIT_INTERVALS={input.probes} "
        "TARGET_INTERVALS={input.hsTarget} "
        "PER_TARGET_COVERAGE={output.target} "
        "COVERAGE_CAP=null"
        ">& {log} "


rule picard_InsertSizeMetrics:
    input:
        bam=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/recalibrated/{sample}.dedup.recal.bam",
        ),
    output:
        metrics=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/recalibrated/{sample}.dedup.recal.is.txt",
        ),
        histogram=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/recalibrated/{sample}.dedup.recal.is.pdf",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/picard.yaml"
        )
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=1),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/picard/InsertsizeMetrics/{sample}.insertsize.log",
        ),
    benchmark:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "benchmarks/picard/IsMetrics/{sample}.txt",
        )
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "picard {params.custom} CollectInsertSizeMetrics "
        "INPUT={input.bam} "
        "OUTPUT={output.metrics} "
        "HISTOGRAM_FILE={output.histogram} "
        ">& {log} "


rule picard_gc_bias:
    input:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "reads/recalibrated/{sample}.dedup.recal.bam",
        ),
    output:
        chart=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "qc/picard/{sample}_gc_bias_metrics.pdf",
        ),
        summary=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "qc/picard/{sample}_summary_metrics.txt",
        ),
        out=resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "qc/picard/{sample}_gc_bias_metrics.txt",
        ),
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=1),
        genome=config.get("resources").get("reference"),
        param=config.get("params").get("picard_gc").get("arguments"),
        tmp_dir=config.get("paths").get("tmp_dir"),
    log:
        resolve_results_filepath(
            config.get("paths").get("results_dir"),
            "logs/picard/CollectGcBiasMetrics/{sample}.gcbias.log",
        ),
    conda:
        resolve_single_filepath(
            config.get("paths").get("workdir"), "workflow/envs/picard.yaml"
        )
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    resources:
        tmpdir=config.get("paths").get("tmp_dir"),
    shell:
        "picard "
        "CollectGcBiasMetrics "
        "{params.custom} "
        "I={input} "
        "R={params.genome} "
        "CHART={output.chart} "
        "S={output.summary} "
        "O={output.out} "
        ">& {log} "
