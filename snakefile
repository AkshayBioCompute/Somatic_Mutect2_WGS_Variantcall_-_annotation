# Snakefile for Somatic Variant Calling using Mutect2

# Define paths
REF_DIR = "/home/akshay/Akshay/Somatic/supporting_files/hg38"
READS_DIR = "/home/akshay/Akshay/Somatic/reads"
ALIGNED_DIR = "/home/akshay/Akshay/Somatic/somatic_mutect2/aligned"
RESULTS_DIR = "/home/akshay/Akshay/Somatic/somatic_mutect2/results"
MUTECT2_FILES_DIR = "/home/akshay/Akshay/Somatic/supporting_files/mutect2_supporting_files"
REF_FASTA = f"{REF_DIR}/hg38.fa"
KNOWN_SITES = f"{REF_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf"
SAMPLES = ["HG008-N-D", "HG008-T"]

# Step 1: QC - FastQC
rule fastqc:
    input:
        r1=f"{READS_DIR}/{sample}.R1.fastq.gz",
        r2=f"{READS_DIR}/{sample}.R2.fastq.gz"
    output:
        r1_out=f"{RESULTS_DIR}/{sample}_R1_fastqc.html",
        r2_out=f"{RESULTS_DIR}/{sample}_R2_fastqc.html"
    shell:
        "fastqc {input.r1} -o {RESULTS_DIR} && fastqc {input.r2} -o {RESULTS_DIR}"

# Step 2: Map to reference using BWA-MEM
rule bwa_mem:
    input:
        r1=f"{READS_DIR}/{sample}.R1.fastq.gz",
        r2=f"{READS_DIR}/{sample}.R2.fastq.gz",
        ref=REF_FASTA
    output:
        f"{ALIGNED_DIR}/{sample}.paired.sam"
    params:
        rg="@RG\\tID:{wildcards.sample}\\tPL:ILLUMINA\\tSM:{wildcards.sample}"
    shell:
        "bwa index {input.ref} && "
        "bwa mem -t 4 -R '{params.rg}' {input.ref} {input.r1} {input.r2} > {output}"

# Step 3: Mark Duplicates
rule mark_duplicates:
    input:
        sam=f"{ALIGNED_DIR}/{sample}.paired.sam"
    output:
        bam=f"{ALIGNED_DIR}/{sample}_sorted_dedup_reads.bam"
    shell:
        "/home/akshay/Akshay/Somatic/tools/gatk/gatk MarkDuplicatesSpark -I {input.sam} -O {output}"

# Step 4: Base Quality Recalibration
rule base_recalibration:
    input:
        bam=f"{ALIGNED_DIR}/{sample}_sorted_dedup_reads.bam",
        ref=REF_FASTA,
        known_sites=KNOWN_SITES
    output:
        recal_table=f"{ALIGNED_DIR}/{sample}_recal_data.table"
    shell:
        "/home/akshay/Akshay/Somatic/tools/gatk/gatk BaseRecalibrator -I {input.bam} -R {input.ref} --known-sites {input.known_sites} -O {output}"

rule apply_bqsr:
    input:
        bam=f"{ALIGNED_DIR}/{sample}_sorted_dedup_reads.bam",
        ref=REF_FASTA,
        recal_table=f"{ALIGNED_DIR}/{sample}_recal_data.table"
    output:
        bqsr_bam=f"{ALIGNED_DIR}/{sample}_sorted_dedup_bqsr_reads.bam"
    shell:
        "/home/akshay/Akshay/Somatic/tools/gatk/gatk ApplyBQSR -I {input.bam} -R {input.ref} --bqsr-recal-file {input.recal_table} -O {output}"

# Step 5: Collect Metrics
rule collect_metrics:
    input:
        bam=f"{ALIGNED_DIR}/{sample}_sorted_dedup_bqsr_reads.bam",
        ref=REF_FASTA
    output:
        alignment_metrics=f"{ALIGNED_DIR}/{sample}_alignment_metrics.txt",
        insert_size_metrics=f"{ALIGNED_DIR}/{sample}_insert_size_metrics.txt",
        insert_size_histogram=f"{ALIGNED_DIR}/{sample}_insert_size_histogram.pdf"
    shell:
        "/home/akshay/Akshay/Somatic/tools/gatk/gatk CollectAlignmentSummaryMetrics R={input.ref} I={input.bam} O={output.alignment_metrics} && "
        "/home/akshay/Akshay/Somatic/tools/gatk/gatk CollectInsertSizeMetrics INPUT={input.bam} OUTPUT={output.insert_size_metrics} HISTOGRAM_FILE={output.insert_size_histogram}"

# Step 6: Call Variants - Mutect2
rule call_variants:
    input:
        normal_bam=f"{ALIGNED_DIR}/HG008-N-D_sorted_dedup_bqsr_reads.bam",
        tumor_bam=f"{ALIGNED_DIR}/HG008-T_sorted_dedup_bqsr_reads.bam",
        ref=REF_FASTA,
        gnomad_vcf=f"{MUTECT2_FILES_DIR}/af-only-gnomad.hg38.vcf.gz",
        pon_vcf=f"{MUTECT2_FILES_DIR}/1000g_pon.hg38.vcf.gz",
        intervals=f"{MUTECT2_FILES_DIR}/exome_calling_regions.v1.1.interval_list"
    output:
        vcf=f"{RESULTS_DIR}/somatic_variants.vcf"
    shell:
        "/home/akshay/Akshay/Somatic/tools/gatk/gatk Mutect2 -R {input.ref} -I {input.tumor_bam} -I {input.normal_bam} "
        "--normal-sample HG008-N-D --tumor-sample HG008-T --germline-resource {input.gnomad_vcf} --panel-of-normals {input.pon_vcf} "
        "-L {input.intervals} -O {output.vcf}"

# Step 7: Estimate Cross-Sample Contamination
rule estimate_contamination:
    input:
        vcf=f"{RESULTS_DIR}/somatic_variants.vcf",
        ref=REF_FASTA
    output:
        contamination_report=f"{RESULTS_DIR}/contamination_report.txt"
    shell:
        "/home/akshay/Akshay/Somatic/tools/gatk/gatk CalculateContamination -I {input.vcf} -R {input.ref} -O {output.contamination_report}"

# Step 8: Estimate Read Orientation Artifacts
rule orientation_artifacts:
    input:
        vcf=f"{RESULTS_DIR}/somatic_variants.vcf",
        ref=REF_FASTA
    output:
        artifacts_report=f"{RESULTS_DIR}/orientation_artifacts_report.txt"
    shell:
        "/home/akshay/Akshay/Somatic/tools/gatk/gatk GetPileupSummaries -I {input.vcf} -R {input.ref} -O {output.artifacts_report}"

# Step 9: Filter Variants
rule filter_variants:
    input:
        vcf=f"{RESULTS_DIR}/somatic_variants.vcf",
        ref=REF_FASTA
    output:
        filtered_vcf=f"{RESULTS_DIR}/filtered_somatic_variants.vcf"
    shell:
        "/home/akshay/Akshay/Somatic/tools/gatk/gatk VariantFiltration -R {input.ref} -V {input.vcf} --filter-expression 'QD < 2.0 || FS > 60.0' --filter-name 'filter1' -O {output.filtered_vcf}"

# Step 10: Annotate Variants
rule annotate_variants:
    input:
        vcf=f"{RESULTS_DIR}/filtered_somatic_variants.vcf",
        ref=REF_FASTA
    output:
        annotated_vcf=f"{RESULTS_DIR}/annotated_variants.vcf"
    shell:
        "/home/akshay/Akshay/Somatic/tools/snpEff/snpEff ann -v GRCh38.99 {input.vcf} > {output.annotated_vcf}"

