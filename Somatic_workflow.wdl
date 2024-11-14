version 1.0

workflow somatic_variant_calling {
  
  input {
    String ref_fasta
    String known_sites
    Array[String] samples
    String gnomad_vcf
    String pon_vcf
    String intervals
  }

  # Step 1: QC - FastQC
  call fastqc

  # Step 2: Map to reference using BWA-MEM
  call bwa_mem

  # Step 3: Mark Duplicates
  call mark_duplicates

  # Step 4: Base Quality Recalibration
  call base_recalibration

  # Step 5: Apply BQSR
  call apply_bqsr

  # Step 6: Collect Metrics
  call collect_metrics

  # Step 7: Call Variants - Mutect2
  call call_variants

  # Step 8: Estimate Cross-Sample Contamination
  call estimate_contamination

  # Step 9: Estimate Read Orientation Artifacts
  call orientation_artifacts

  # Step 10: Filter Variants
  call filter_variants

  # Step 11: Annotate Variants
  call annotate_variants
}

task fastqc {
  input {
    File r1
    File r2
  }

  command {
    fastqc ${r1} -o . && fastqc ${r2} -o .
  }

  output {
    File r1_out = "${r1}.html"
    File r2_out = "${r2}.html"
  }
}

task bwa_mem {
  input {
    File r1
    File r2
    File ref_fasta
  }

  command {
    bwa mem -t 4 -R "@RG\\tID:${sample}\\tPL:ILLUMINA\\tSM:${sample}" ${ref_fasta} ${r1} ${r2} > ${sample}.paired.sam
  }

  output {
    File sam_out = "${sample}.paired.sam"
  }
}

task mark_duplicates {
  input {
    File sam
  }

  command {
    gatk MarkDuplicatesSpark -I ${sam} -O ${sam}.dedup.bam
  }

  output {
    File bam_out = "${sam}.dedup.bam"
  }
}

task base_recalibration {
  input {
    File bam
    File ref_fasta
    File known_sites
  }

  command {
    gatk BaseRecalibrator -I ${bam} -R ${ref_fasta} --known-sites ${known_sites} -O ${bam}.recal.table
  }

  output {
    File recal_table = "${bam}.recal.table"
  }
}

task apply_bqsr {
  input {
    File bam
    File ref_fasta
    File recal_table
  }

  command {
    gatk ApplyBQSR -I ${bam} -R ${ref_fasta} --bqsr-recal-file ${recal_table} -O ${bam}.bqsr.bam
  }

  output {
    File bqsr_bam = "${bam}.bqsr.bam"
  }
}

task collect_metrics {
  input {
    File bam
    File ref_fasta
  }

  command {
    gatk CollectAlignmentSummaryMetrics -R ${ref_fasta} -I ${bam} -O ${bam}.alignment_metrics.txt
    gatk CollectInsertSizeMetrics -I ${bam} -O ${bam}.insert_size_metrics.txt -H ${bam}.insert_size_histogram.pdf
  }

  output {
    File alignment_metrics = "${bam}.alignment_metrics.txt"
    File insert_size_metrics = "${bam}.insert_size_metrics.txt"
    File insert_size_histogram = "${bam}.insert_size_histogram.pdf"
  }
}

task call_variants {
  input {
    File normal_bam
    File tumor_bam
    File ref_fasta
    File gnomad_vcf
    File pon_vcf
    File intervals
  }

  command {
    gatk Mutect2 -R ${ref_fasta} -I ${tumor_bam} -I ${normal_bam} --normal-sample HG008-N-D --tumor-sample HG008-T \
      --germline-resource ${gnomad_vcf} --panel-of-normals ${pon_vcf} -L ${intervals} -O somatic_variants.vcf
  }

  output {
    File vcf = "somatic_variants.vcf"
  }
}

task estimate_contamination {
  input {
    File vcf
    File ref_fasta
  }

  command {
    gatk CalculateContamination -I ${vcf} -R ${ref_fasta} -O contamination_report.txt
  }

  output {
    File contamination_report = "contamination_report.txt"
  }
}

task orientation_artifacts {
  input {
    File vcf
    File ref_fasta
  }

  command {
    gatk GetPileupSummaries -I ${vcf} -R ${ref_fasta} -O orientation_artifacts_report.txt
  }

  output {
    File artifacts_report = "orientation_artifacts_report.txt"
  }
}

task filter_variants {
  input {
    File vcf
    File ref_fasta
  }

  command {
    gatk VariantFiltration -R ${ref_fasta} -V ${vcf} --filter-expression 'QD < 2.0 || FS > 60.0' --filter-name 'filter1' -O filtered_somatic_variants.vcf
  }

  output {
    File filtered_vcf = "filtered_somatic_variants.vcf"
  }
}

task annotate_variants {
  input {
    File vcf
    File ref_fasta
  }

  command {
    snpEff ann -v GRCh38.99 ${vcf} > annotated_variants.vcf
  }

  output {
    File annotated_vcf = "annotated_variants.vcf"
  }
}

