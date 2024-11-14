#!/usr/bin/env nextflow

// Define paths
params.ref_dir = "/home/akshay/Akshay/Somatic/supporting_files/hg38"
params.reads_dir = "/home/akshay/Akshay/Somatic/reads"
params.aligned_dir = "/home/akshay/Akshay/Somatic/somatic_mutect2/aligned"
params.results_dir = "/home/akshay/Akshay/Somatic/somatic_mutect2/results"
params.mutect2_files_dir = "/home/akshay/Akshay/Somatic/supporting_files/mutect2_supporting_files"
params.ref_fasta = "${params.ref_dir}/hg38.fa"
params.known_sites = "${params.ref_dir}/Homo_sapiens_assembly38.dbsnp138.vcf"
params.samples = ["HG008-N-D", "HG008-T"]

process fastqc {
    input:
    tuple val(sample)

    output:
    tuple path("${params.results_dir}/${sample}_R1_fastqc.html"), path("${params.results_dir}/${sample}_R2_fastqc.html")

    script:
    """
    fastqc ${params.reads_dir}/${sample}.R1.fastq.gz -o ${params.results_dir}
    fastqc ${params.reads_dir}/${sample}.R2.fastq.gz -o ${params.results_dir}
    """
}

process bwa_mem {
    input:
    tuple val(sample)

    output:
    path("${params.aligned_dir}/${sample}.paired.sam")

    script:
    """
    bwa index ${params.ref_fasta}
    bwa mem -t 4 -R "@RG\\tID:${sample}\\tPL:ILLUMINA\\tSM:${sample}" ${params.ref_fasta} ${params.reads_dir}/${sample}.R1.fastq.gz ${params.reads_dir}/${sample}.R2.fastq.gz > ${params.aligned_dir}/${sample}.paired.sam
    """
}

process mark_duplicates {
    input:
    path("${params.aligned_dir}/${sample}.paired.sam")

    output:
    path("${params.aligned_dir}/${sample}_sorted_dedup_reads.bam")

    script:
    """
    /home/akshay/Akshay/Somatic/tools/gatk/gatk MarkDuplicatesSpark -I ${params.aligned_dir}/${sample}.paired.sam -O ${params.aligned_dir}/${sample}_sorted_dedup_reads.bam
    """
}

process base_recalibration {
    input:
    path("${params.aligned_dir}/${sample}_sorted_dedup_reads.bam"),
    path("${params.ref_fasta}"),
    path("${params.known_sites}")

    output:
    path("${params.aligned_dir}/${sample}_recal_data.table")

    script:
    """
    /home/akshay/Akshay/Somatic/tools/gatk/gatk BaseRecalibrator -I ${params.aligned_dir}/${sample}_sorted_dedup_reads.bam -R ${params.ref_fasta} --known-sites ${params.known_sites} -O ${params.aligned_dir}/${sample}_recal_data.table
    """
}

process apply_bqsr {
    input:
    path("${params.aligned_dir}/${sample}_sorted_dedup_reads.bam"),
    path("${params.ref_fasta}"),
    path("${params.aligned_dir}/${sample}_recal_data.table")

    output:
    path("${params.aligned_dir}/${sample}_sorted_dedup_bqsr_reads.bam")

    script:
    """
    /home/akshay/Akshay/Somatic/tools/gatk/gatk ApplyBQSR -I ${params.aligned_dir}/${sample}_sorted_dedup_reads.bam -R ${params.ref_fasta} --bqsr-recal-file ${params.aligned_dir}/${sample}_recal_data.table -O ${params.aligned_dir}/${sample}_sorted_dedup_bqsr_reads.bam
    """
}

process collect_metrics {
    input:
    path("${params.aligned_dir}/${sample}_sorted_dedup_bqsr_reads.bam"),
    path("${params.ref_fasta}")

    output:
    path("${params.aligned_dir}/${sample}_alignment_metrics.txt"),
    path("${params.aligned_dir}/${sample}_insert_size_metrics.txt"),
    path("${params.aligned_dir}/${sample}_insert_size_histogram.pdf")

    script:
    """
    /home/akshay/Akshay/Somatic/tools/gatk/gatk CollectAlignmentSummaryMetrics R=${params.ref_fasta} I=${params.aligned_dir}/${sample}_sorted_dedup_bqsr_reads.bam O=${params.aligned_dir}/${sample}_alignment_metrics.txt
    /home/akshay/Akshay/Somatic/tools/gatk/gatk CollectInsertSizeMetrics INPUT=${params.aligned_dir}/${sample}_sorted_dedup_bqsr_reads.bam OUTPUT=${params.aligned_dir}/${sample}_insert_size_metrics.txt HISTOGRAM_FILE=${params.aligned_dir}/${sample}_insert_size_histogram.pdf
    """
}

process call_variants {
    input:
    path("${params.aligned_dir}/HG008-N-D_sorted_dedup_bqsr_reads.bam"),
    path("${params.aligned_dir}/HG008-T_sorted_dedup_bqsr_reads.bam"),
    path("${params.ref_fasta}"),
    path("${params.mutect2_files_dir}/af-only-gnomad.hg38.vcf.gz"),
    path("${params.mutect2_files_dir}/1000g_pon.hg38.vcf.gz"),
    path("${params.mutect2_files_dir}/exome_calling_regions.v1.1.interval_list")

    output:
    path("${params.results_dir}/somatic_variants.vcf")

    script:
    """
    /home/akshay/Akshay/Somatic/tools/gatk/gatk Mutect2 -R ${params.ref_fasta} -I ${params.aligned_dir}/HG008-T_sorted_dedup_bqsr_reads.bam -I ${params.aligned_dir}/HG008-N-D_sorted_dedup_bqsr_reads.bam --normal-sample HG008-N-D --tumor-sample HG008-T --germline-resource ${params.mutect2_files_dir}/af-only-gnomad.hg38.vcf.gz --panel-of-normals ${params.mutect2_files_dir}/1000g_pon.hg38.vcf.gz -L ${params.mutect2_files_dir}/exome_calling_regions.v1.1.interval_list -O ${params.results_dir}/somatic_variants.vcf
    """
}

process estimate_contamination {
    input:
    path("${params.results_dir}/somatic_variants.vcf"),
    path("${params.ref_fasta}")

    output:
    path("${params.results_dir}/contamination_report.txt")

    script:
    """
    /home/akshay/Akshay/Somatic/tools/gatk/gatk CalculateContamination -I ${params.results_dir}/somatic_variants.vcf -R ${params.ref_fasta} -O ${params.results_dir}/contamination_report.txt
    """
}

process orientation_artifacts {
    input:
    path("${params.results_dir}/somatic_variants.vcf"),
    path("${params.ref_fasta}")

    output:
    path("${params.results_dir}/orientation_artifacts_report.txt")

    script:
    """
    /home/akshay/Akshay/Somatic/tools/gatk/gatk GetPileupSummaries -I ${params.results_dir}/somatic_variants.vcf -R ${params.ref_fasta} -O ${params.results_dir}/orientation_artifacts_report.txt
    """
}

process filter_variants {
    input:
    path("${params.results_dir}/somatic_variants.vcf"),
    path("${params.ref_fasta}")

    output:
    path("${params.results_dir}/filtered_somatic_variants.vcf")

    script:
    """
    /home/akshay/Akshay/Somatic/tools/gatk/gatk VariantFiltration -R ${params.ref_fasta} -V ${params.results_dir}/somatic_variants.vcf --filter-expression 'QD < 2.0 || FS > 60.0' --filter-name 'filter1' -O ${params.results_dir}/filtered_somatic_variants.vcf
    """
}

process annotate_variants {
    input:
    path("${params.results_dir}/filtered_somatic_variants.vcf"),
    path("${params.ref_fasta}")

    output:
    path("${params.results_dir}/annotated_variants.vcf")

    script:
    """
    /home/akshay/Akshay/Somatic/tools/snpEff/snpEff ann -v GRCh38.99 ${params.results_dir}/filtered_somatic_variants.vcf > ${params.results_dir}/annotated_variants.vcf
    """
}

workflow {
    params.samples.each { sample ->
        fastqc(sample)
        bwa_mem(sample)
        mark_duplicates(sample)
        base_recalibration(sample)
        apply_bqsr(sample)
        collect_metrics(sample)
    }

    call_variants()
    estimate_contamination()
    orientation_artifacts()
    filter_variants()
    annotate_variants()
}

