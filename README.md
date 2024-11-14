Here’s a sample `README.md` for your project, covering the Snakefile, WDL, Nextflow, requirements, and details about your `src` folder.

---

# Genomic Analysis : Somatic_Mutect2_Variantcall_-_annotation

## Overview

This repository provides a genomic analysis pipeline that includes the necessary tools and scripts to process sequencing data, from raw reads through alignment, variant calling, and annotation. The pipeline is designed to support multiple workflow management systems (WMS), including SnakeMake, WDL, and Nextflow. It also includes a set of Python scripts for additional utility functions and quality control.

---

## Table of Contents

- [Workflow Management Systems (WMS)](#workflow-management-systems-wms)
  - [SnakeMake](#snakefile)
  - [WDL](#wdl)
  - [Nextflow](#nextflow)
- [Requirements](#requirements)
- [src Folder](#src-folder)
- [Getting Started](#getting-started)
  - [Installation](#installation)
  - [Running the Pipeline](#running-the-pipeline)
- [Contributing](#contributing)
- [License](#license)

---

## Workflow Management Systems (WMS)

### SnakeMake

The SnakeMake workflow (`Snakefile`) is designed to orchestrate the entire pipeline. SnakeMake is a powerful tool to manage computational workflows. The `Snakefile` defines the rules and dependencies between different pipeline steps.

### WDL

WDL (Workflow Description Language) scripts are provided to run the same pipeline on cloud platforms such as Google Cloud and AWS. The WDL files can be executed with platforms such as [Cromwell](https://cromwell.readthedocs.io/).

### Nextflow

Nextflow scripts are included to run the pipeline in a highly parallelized and containerized manner, making it easy to run the pipeline on a cluster or cloud environment.

---

## Requirements

The following Python packages and software dependencies are required to run this pipeline:

- **Python packages**:  
  `requests`, `vcf`, `fastqc`, and any other necessary dependencies are listed in the `requirements.txt`.

- **Bioinformatics tools**:
  - BWA
  - GATK
  - SnpEff
  - FastQC
  - Picard (for MarkDuplicates)

### `requirements.txt`
This file contains the Python dependencies required for utility scripts.

```txt
requests==2.25.1
vcf==0.8.8
fastqc==0.11.9
```

You can install these packages with:

```bash
pip install -r requirements.txt
```

---

## src Folder

The `src` folder contains all the necessary scripts for aligning reads, processing BAM files, variant calling, and annotating variants. Key files in the `src` folder include:

- **Alignment Scripts**:
  - `bwa_mem.sh`: Aligns paired-end reads using BWA MEM.
  - `bwa_index.sh`: Creates BWA index for reference genome.

- **GATK Scripts**:
  - `mark_duplicates.sh`: Marks duplicate reads in BAM files.
  - `base_recalibration.sh`: Performs base quality score recalibration.
  - `apply_bqsr.sh`: Applies base quality score recalibration.
  - `call_variants.sh`: Calls variants using Mutect2.
  - `filter_variants.sh`: Filters variants based on quality.
  - `annotate_variants.sh`: Annotates variants using SnpEff.

- **Utility Scripts**:
  - `download_data.py`: Downloads reference genome and other necessary files.
  - `process_vcf.py`: Filters and processes VCF files.
  - `qc_tools.py`: Parses and processes FastQC reports.

---

## Getting Started

### Installation

To run the pipeline, first, install the necessary tools and dependencies.

1. **Clone the repository**:
   ```bash
   git clone <repository_url>
   cd <repository_directory>
   ```

2. **Install Python dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

3. **Install bioinformatics tools**:
   Ensure that you have installed the following tools on your system:
   - BWA
   - GATK
   - SnpEff
   - Picard (MarkDuplicates)
   - FastQC

4. **Ensure you have the correct reference genome**:
   The pipeline requires a reference genome (e.g., GRCh38 for human). Make sure you have it downloaded and indexed.

---

### Running the Pipeline

The pipeline can be run using one of the following workflow management systems:

#### Using SnakeMake:

1. **Run SnakeMake**:
   ```bash
   snakemake -j 4  # Adjust the number of jobs as needed
   ```

#### Using WDL (Cromwell):

1. **Execute with Cromwell**:
   - Install [Cromwell](https://cromwell.readthedocs.io/).
   - Run the WDL workflow:
     ```bash
     java -jar cromwell.jar run <workflow.wdl>
     ```

#### Using Nextflow:

1. **Run Nextflow**:
   - Install [Nextflow](https://www.nextflow.io/).
   - Execute the pipeline:
     ```bash
     nextflow run <pipeline.nf>
     ```

---

## Contributing

We welcome contributions! If you'd like to improve the pipeline or add new features, please:

1. Fork the repository
2. Create a new branch
3. Make your changes
4. Submit a pull request with a description of your changes

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

This `README.md` provides a comprehensive overview of how to set up, install, and run the genomic analysis pipeline with various workflow management systems. You can adjust this document to suit specific details or configurations that may differ in your setup.
