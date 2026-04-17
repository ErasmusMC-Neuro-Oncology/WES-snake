# README.md

## Overview

This pipeline performs **tumor-only somatic variant calling and copy number analysis on FFPE tumor material** using a Snakemake workflow. It is designed for targeted or exome sequencing data where no matched normal sample is available.

The workflow combines:

* **nf-core/Sarek** for alignment and mutation calling with Mutect2
* **CopywriteR / QDNAseq** for copy number profiling
* **PureCN** for tumor-only purity/ploidy estimation and somatic vs germline classification
* Summary reporting and QC plots

The main goal is to generate a high-confidence list of somatic variants and copy number alterations from FFPE tumor sequencing data.

---

## Workflow Summary

### 1. Sample Preparation

The pipeline first creates a `samplesheet.csv` from input metadata files. This script has to be edited by the user

### 2. Preprocessing and mutation calling with Sarek

For each patient:

* Reads are aligned to the reference genome
* Duplicate marking and recalibration are performed
* Variants are called using **Mutect2**
* Variants are annotated
* Additional annotation with:

  * Panel of Normals (PON)
  * COSMIC

### 3. Copy Number Analysis

Using tumor BAM files:

* Off-target read depth CNA calling with **CopywriteR**
* Segmentation with **QDNAseq**
* Copy number calling
* Purity/ploidy fitting using **ACE**
* Homologous recombination deficiency metrics using **CNH**

### 4. Tumor-Only Interpretation (PureCN)

PureCN uses:

* Variant allele frequencies
  n- Coverage
* Segmentation data

To estimate:

* Tumor purity
* Tumor ploidy
* Somatic vs germline probability
* Cancer cell fraction (CCF)
* Tumor mutational burden (TMB)
* Mutational signatures

### 5. Reporting

Merged summary tables and QC plots are generated.

---

## Software Requirements

Tested with:

* **Snakemake 9.17.3**
* **nf-core/Sarek > 3.8.1**
* Conda / Mamba
* Samtools
* bcftools
* R

---

## Running the Pipeline

### Option 1 — SLURM Cluster

```bash
snakemake --profile slurm/
```

### Option 2 — Local Workstation

```bash
snakemake --cores 1 --use-conda
```

---

## Main Output Files

### Final Variant Calls

```text
PureCN/{binsize}/{patient}/{patient}_tumor1_variants.csv
```

Primary final mutation output flagged by somatic status.

### Final Copy Number Outputs

```text
QDNAseq/{binsize}/{patient}/data/QDNAseq_Segments.txt
QDNAseq/{binsize}/{patient}/data/QDNAseq_calls.txt
```

Primary CNA outputs.

### Additional Outputs

#### Sample Overview Table

```text
sampledata/SampleData_WES.txt
```

A summary table containing per-sample overview metrics such as sequencing depth, purity estimates, CNA statistics, mutational burden, and other aggregated QC/output metrics.

```text
PureCN/{binsize}/{patient}/{patient}_tumor1.csv
PureCN/{binsize}/{patient}/{patient}_tumor1_mutation_burden.csv
PureCN/{binsize}/{patient}/{patient}_tumor1_signatures.csv
sarek/{patient}/annotation/mutect2/{patient}_tumor1/{patient}_tumor1.annotated.vcf.gz
```

---

## Directory Structure

```text
output/
├── sarek/
├── PureCN/
├── QDNAseq/
├── ACE/
├── CNH/
├── sampledata/
└── plots/
```

---

## Notes on FFPE Samples

Because FFPE material can contain artifacts:

* Mutect2 filtering is applied
* PureCN assists somatic vs germline classification
* Copy number calls should be interpreted with tumor purity estimates

**By default, IDH1/2 hotspot mutations are given a high (0.999) prior somatic probability. This is done in the scripts/PureCN.R script.**
