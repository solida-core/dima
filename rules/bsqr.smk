def get_known_sites(known_sites=['dbsnp','mills','ph1_indel']):
    known_variants = config.get("known_variants")
    ks = []
    if len(known_sites) == 0:
        known_sites = known_variants.keys()
    for k, v in known_variants.items():
        if k in known_sites:
            ks.append("--known-sites {} ".format(resolve_single_filepath(
                *references_abs_path(), v)[0]))
    return "".join(ks)


rule gatk_BQSR_data_processing:
    input:
        cram="reads/dedup/{sample}.dedup.cram"
    output:
        "reads/recalibrated/{sample}.recalibrate.grp"
    conda:
       "../envs/gatk.yaml"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        known_sites=get_known_sites(config.get("rules").get("gatk_BQSR").get("known_sites"))
    log:
        "logs/gatk/BaseRecalibrator/{sample}_BQSR_data_processing_info.log"
    benchmark:
        "benchmarks/gatk/BaseRecalibrator/{sample}_BaseRecalibrator_data_processing_info.txt"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk BaseRecalibrator --java-options {params.custom} "
        "-R {params.genome} "
        "{params.known_sites} "
#        "--spark-runner LOCAL "
#        "--spark-master local[{threads}]  "
        "-I {input.cram} "
        "-O {output} "
        ">& {log}"


rule gatk_ApplyBQSR:
    input:
        cram="reads/dedup/{sample}.dedup.cram",
        bqsr="reads/recalibrated/{sample}.recalibrate.grp"
    output:
        cram=protected("reads/recalibrated/{sample}.dedup.recal.cram"),
        crai="reads/recalibrated/{sample}.dedup.recal.crai"
    conda:
        "../envs/gatk.yaml"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta"))
    log:
        "logs/gatk/ApplyBQSR/{sample}.post_recalibrate_info.log"
    benchmark:
        "benchmarks/gatk/ApplyBQSR/{sample}.post_recalibrate_info.txt"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk ApplyBQSR --java-options {params.custom} "
        "-R {params.genome} "
#        "--spark-runner LOCAL "
#        "--spark-master local[{threads}]  "
        "-I {input.cram} "
        "--bqsr-recal-file {input.bqsr} "
        "-O {output.cram} "
        ">& {log}"


rule gatk_BQSR_quality_control:
    input:
        cram="reads/recalibrated/{sample}.dedup.recal.cram",
        pre="reads/recalibrated/{sample}.recalibrate.grp"
    output:
        post="reads/recalibrated/{sample}.post.recalibrate.grp",
        plot="reads/recalibrated/{sample}.recalibration_plots.pdf"
    conda:
        "../envs/gatk.yaml"
    params:
        custom=java_params(tmp_dir=config.get("tmp_dir"), multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        known_sites=get_known_sites(config.get("rules").get("gatk_BQSR").get("known_sites"))
    log:
        b="logs/gatk/BaseRecalibrator/{sample}_BQSR_quality_control_info.log",
        a="logs/gatk/AnalyzeCovariates/{sample}_BQSR_quality_control_cov_info.log"
    benchmark:
        "benchmarks/gatk/BaseRecalibrator/{sample}_BQSR_quality_control_info.txt"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk BaseRecalibrator --java-options {params.custom} "
        "-R {params.genome} "
        "{params.known_sites} "
#        "--spark-runner LOCAL "
#        "--spark-master local[{threads}]  "
        "-I {input.cram} "
        "-O {output.post} "
        ">& {log.b}; "
        "gatk AnalyzeCovariates --java-options {params.custom} "
        "-before {input.pre} -after {output.post} -plots {output.plot} "
        ">& {log.a}"
