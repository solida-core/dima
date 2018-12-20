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
        bam=dima_path+"reads/dedup/{sample}.dedup.bam"
    output:
        dima_path+"reads/recalibrated/{sample}.recalibrate.grp"
    conda:
       "../envs/gatk.yaml"
    params:
        custom=java_params(tmp_dir=config.get("paths").get("to_tmp"),
                           multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        known_sites=get_known_sites(config.get("rules").get("gatk_BQSR").get("known_sites"))
    log:
        dima_path+"logs/gatk/BaseRecalibrator/{sample}_BQSR_data_processing_info.log"
    benchmark:
        dima_path+"benchmarks/gatk/BaseRecalibrator/{sample}_BaseRecalibrator_data_processing_info.txt"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk BaseRecalibrator --java-options {params.custom} "
        "-R {params.genome} "
        "{params.known_sites} "
#        "--spark-runner LOCAL "
#        "--spark-master local[{threads}]  "
        "-I {input.bam} "
        "-O {output} "
        ">& {log}"


rule gatk_ApplyBQSR:
    input:
        bam=dima_path+"reads/dedup/{sample}.dedup.bam",
        bqsr=dima_path+"reads/recalibrated/{sample}.recalibrate.grp"
    output:
        bam=protected(dima_path+"reads/recalibrated/{sample}.dedup.recal.bam"),

    conda:
        "../envs/gatk.yaml"
    params:
        custom=java_params(tmp_dir=config.get("paths").get("to_tmp"),
                           multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta"))
    log:
        dima_path+"logs/gatk/ApplyBQSR/{sample}.post_recalibrate_info.log"
    benchmark:
        dima_path+"benchmarks/gatk/ApplyBQSR/{sample}.post_recalibrate_info.txt"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk ApplyBQSR --java-options {params.custom} "
        "-R {params.genome} "
#        "--spark-runner LOCAL "
#        "--spark-master local[{threads}]  "
        "-I {input.bam} "
        "--bqsr-recal-file {input.bqsr} "
        "-O {output.bam} "
        ">& {log}"


rule gatk_BQSR_quality_control:
    input:
        bam=dima_path+"reads/recalibrated/{sample}.dedup.recal.bam",
        pre=dima_path+"reads/recalibrated/{sample}.recalibrate.grp"
    output:
        post=dima_path+"reads/recalibrated/{sample}.post.recalibrate.grp",
        plot=dima_path+"reads/recalibrated/{sample}.recalibration_plots.pdf"
    conda:
        "../envs/gatk.yaml"
    params:
        custom=java_params(tmp_dir=config.get("paths").get("to_tmp"),
                           multiply_by=5),
        genome=resolve_single_filepath(*references_abs_path(), config.get("genome_fasta")),
        known_sites=get_known_sites(config.get("rules").get("gatk_BQSR").get("known_sites"))
    log:
        b=dima_path+"logs/gatk/BaseRecalibrator/{sample}_BQSR_quality_control_info.log",
        a=dima_path+"logs/gatk/AnalyzeCovariates/{sample}_BQSR_quality_control_cov_info.log"
    benchmark:
        dima_path+"benchmarks/gatk/BaseRecalibrator/{sample}_BQSR_quality_control_info.txt"
    threads: conservative_cpu_count(reserve_cores=2, max_cores=99)
    shell:
        "gatk BaseRecalibrator --java-options {params.custom} "
        "-R {params.genome} "
        "{params.known_sites} "
#        "--spark-runner LOCAL "
#        "--spark-master local[{threads}]  "
        "-I {input.bam} "
        "-O {output.post} "
        ">& {log.b}; "
        "gatk AnalyzeCovariates --java-options {params.custom} "
        "-before {input.pre} -after {output.post} -plots {output.plot} "
        ">& {log.a}"
