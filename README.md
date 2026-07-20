# MEL-Targeted-Variant-Analysis

> A reproducible sequencing pipeline for germline and somatic variant analysis in melanoma using targeted gene panels

---

## Overview

This repository contains an end-to-end bioinformatics workflow for targeted panel sequencing analysis of melanoma samples.

The pipeline performs:

- Quality control of raw FASTQ files
- Adapter and quality trimming
- Read alignment to the human reference genome (GRCh38)
- BAM processing and quality assessment
- Germline variant calling
- Somatic variant calling using matched tumour-normal samples
- Functional variant annotation
- Variant interpretation

---

## Workflow

```text
FASTQ
   │
   ▼
Quality Control (FastQC)
   │
   ▼
Read Trimming (fastp / Trimmomatic)
   │
   ▼
Alignment (BWA-MEM)
   │
   ▼
Sorted BAM
   │
   ▼
Duplicate Marking
   │
   ▼
Base Quality Score Recalibration (BQSR)
   │
   ├───────────────┐
   ▼               ▼
Germline         Somatic
HaplotypeCaller  Mutect2
   │               │
   ▼               ▼
     Filtered VCFs
          │
          ▼
  Annotation (Funcotator)
          │
          ▼
 Final Annotated Variants

```
---

## Project Structure

```
MEL-Targeted-Variant-Analysis/
│
├── Data/
│   ├── Reference/
│   ├── Targets/
│
├── Results/
│   ├── QC/
│   ├── BAM/
│   ├── Germline/
│   ├── Somatic/
│   └── Annotation/
│
├── workflow_scripts/
│   ├── 01_GERMLINE_variant_calling_pipeline.sh
│   ├── 02_SOMATIC_variant_calling_pipeline_tumour_vs_normal.sh
│
├── README.md
└── LICENSE
```

## Sequencing Panels

### Germline

| Panel | Description |
|---------|-------------|
| Hereditary Cancer Panel | 58 genes |

### Somatic

| Panel | Description |
|---------|-------------|
| Ion AmpliSeq Melanoma Panel | Targeted melanoma-associated genes |


## Software

| Tool | Purpose |
|------|---------|
| FastQC | Raw read quality assessment |
| fastp / Trimmomatic | Read trimming |
| BWA-MEM | Read alignment |
| Samtools | BAM processing |
| Picard | Duplicate marking |
| GATK | Variant calling |
| Mutect2 | Somatic variant detection |
| HaplotypeCaller | Germline variant detection |
| Funcotator | Variant annotation |

---

## Citation

If you use this workflow in your research, please cite this repository.

---

## License

Released under the MIT License.
