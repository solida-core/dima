
rule gatk_BQSR_data_processing:
    input:
        bam=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/dedup/{sample}.dedup.bam")
    output:
        resolve_results_filepath(config.get("paths").get("results_dir"),"reads/recalibrated/{sample}.recalibrate.grp")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/gatk.yaml")
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=1),
        genome=config.get("resources").get("reference"),
        known_sites=get_known_sites(config.get("params").get("gatk_BQSR").get("known_sites"))
    log:
        resolve_results_filepath(config.get("paths").get("results_dir"),"logs/gatk/BaseRecalibrator/{sample}_BQSR_data_processing_info.log")
    benchmark:
        resolve_results_filepath(config.get("paths").get("results_dir"),"benchmarks/gatk/BaseRecalibrator/{sample}_BaseRecalibrator_data_processing_info.txt")
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk BaseRecalibrator --java-options {params.custom} "
        "-R {params.genome} "
        "{params.known_sites} "
        "-I {input.bam} "
        "-O {output} "
        ">& {log}"


rule gatk_ApplyBQSR:
    input:
        bam=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/dedup/{sample}.dedup.bam"),
        bqsr=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/recalibrated/{sample}.recalibrate.grp")
    output:
        bam=protected(resolve_results_filepath(config.get("paths").get("results_dir"),"reads/recalibrated/{sample}.dedup.recal.bam")),
        bai=protected(resolve_results_filepath(config.get("paths").get("results_dir"),"reads/recalibrated/{sample}.dedup.recal.bai"))
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/gatk.yaml")
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=1),
        genome=config.get("resources").get("reference")
    log:
        resolve_results_filepath(config.get("paths").get("results_dir"),"logs/gatk/ApplyBQSR/{sample}.post_recalibrate_info.log")
    benchmark:
        resolve_results_filepath(config.get("paths").get("results_dir"),"benchmarks/gatk/ApplyBQSR/{sample}.post_recalibrate_info.txt")
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk ApplyBQSR --java-options {params.custom} "
        "-R {params.genome} "
        "-I {input.bam} "
        "--bqsr-recal-file {input.bqsr} "
        "-O {output.bam} "
        ">& {log}"


rule gatk_BQSR_quality_control:
    input:
        bam=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/recalibrated/{sample}.dedup.recal.bam"),
        pre=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/recalibrated/{sample}.recalibrate.grp")
    output:
        post=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/recalibrated/{sample}.post.recalibrate.grp"),
        plot=resolve_results_filepath(config.get("paths").get("results_dir"),"reads/recalibrated/{sample}.recalibration_plots.pdf")
    conda:
        resolve_single_filepath(config.get("paths").get("workdir"),"workflow/envs/gatk.yaml")
    params:
        custom=java_params(tmp_dir=config.get("paths").get("tmp_dir"), multiply_by=1),
        genome=config.get("resources").get("reference"),
        known_sites=get_known_sites(config.get("params").get("gatk_BQSR").get("known_sites"))
    log:
        b=resolve_results_filepath(config.get("paths").get("results_dir"),"logs/gatk/BaseRecalibrator/{sample}_BQSR_quality_control_info.log"),
        a=resolve_results_filepath(config.get("paths").get("results_dir"),"logs/gatk/AnalyzeCovariates/{sample}_BQSR_quality_control_cov_info.log")
    benchmark:
        resolve_results_filepath(config.get("paths").get("results_dir"),"benchmarks/gatk/BaseRecalibrator/{sample}_BQSR_quality_control_info.txt")
    resources:
        tmpdir = config.get("paths").get("tmp_dir")
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk BaseRecalibrator --java-options {params.custom} "
        "-R {params.genome} "
        "{params.known_sites} "
        "-I {input.bam} "
        "-O {output.post} "
        ">& {log.b}; "
        "gatk AnalyzeCovariates --java-options {params.custom} "
        "-before {input.pre} -after {output.post} -plots {output.plot} "
        ">& {log.a}"
