# 🧬 End-to-End Genomic Driver Analysis in Breast Cancer Organoids (PDOs)

## 📌 Project Overview

This repository contains the **complete bioinformatics pipeline** used to analyze Copy Number Variations (CNV) and Somatic Point Mutations (SNV/Indel) for a panel of clinically relevant genes in three Patient-Derived Organoid (PDO) samples (Nr1, Nr2, Nr3).

The primary goal was to generate a **rigorous and biologically coherent OncoPrint** by applying systematic quality filters across three main stages:
1.  **Primary/Secondary Analysis:** QC, Alignment, Variant Calling (GATK).
2.  **Intermediate Processing:** VCF/SEG filtering and summarization (Python).
3.  **Final Visualization:** Biological consistency QC and OncoPrint generation (R).

## Part I: Primary & Secondary Analysis (`scripts/`)

This section details the execution sequence of the core genomic pipeline steps (using SLURM and Singularity images):

| Script | Function | Key Tools |
| :--- | :--- | :--- |
| `01_fastqc.sh` | Initial Quality Control (QC) of raw FASTQ reads. | FastQC |
| `02_trimming.sh` | Adapter and low-quality sequence trimming. | TrimGalore |
| `03_alignment.sh` | Alignment of trimmed reads to the reference genome (GRCh37/hg19). | BWA, Samtools |
| `04_dedup.sh` | Marking of PCR duplicates. | Picard MarkDuplicates |
| `04b_add_readgroup.sh` | Adding Read Group (RG) information, mandatory for GATK tools. | Picard AddOrReplaceReadGroups |
| `05_mutect.sh` | Somatic SNV/Indel detection (tumor-only mode) and initial filtering. | GATK MuTect2, GATK FilterMutectCalls |
| `06_CNV.sh` | CNV segmentation and modeling using GATK CNV. | GATK CNV Suite |

## Part II: Intermediate Processing (`scripts/` - Python)

These Python scripts bridge the gap between the raw genomic output (VCF, SEG) and the structured input required for the final R analysis, using the gene list defined in `resources/genes_interes_2.bed`.

| Script | Function | Output |
| :--- | :--- | :--- |
| `summary_variants2.py` | Filters the annotated VCFs for target genes and summarizes functional impact (HIGH/MODERATE). | `data/processed/variant_summary_genes_2.xlsx` |
| `cnv_gene_3.py` | Processes CNV segments, applies an initial numeric threshold ($\lvert Segment\_Mean \rvert \geq 0.3$), and maps segments to target genes. | `data/processed/CNV_genes_overlap_3.csv` |

## Part III: Final Biological QC and Visualization (`src/07_oncoprint_visualization.R`)

The final R script loads the processed data, implements the most rigorous filtering steps, and generates the final OncoPrint figure.

### A. Strict Numeric CNV Filtering
Only alterations exceeding a robust significance threshold are considered:
* **Criterion:** $\lvert Segment\_Mean \rvert \geq 0.3$.

### B. Biological Consistency Filter (CNV Artifact Removal)

This critical QC step removes alterations with a high probability of being **False Positives** (noise artifacts) by checking consistency against known cancer biology:
* **Exclusion of Deletions in Oncogenes:** Deletions are discarded for classical oncogenes (e.g., **PIK3CA, KRAS, ESR1**), as their *driver* mechanism is activation.
* **Exclusion of Amplifications in Tumor Suppressors:** Amplifications are discarded for tumor suppressors (e.g., **TP53, PTEN**), as their *driver* mechanism is inactivation/loss-of-function.

### C. Driver Variant Rescue and Prioritization

* **PIK3CA Rescue:** The critical *hotspot* mutation **PIK3CA p.E545K** was manually included as it was initially filtered by MuTect2 due to a low Allele Frequency (AF).
* **Prioritization:** In cases where a single gene has both an SNV/Indel and a CNV, the **SNV/Indel is given absolute priority** for visualization.

## Execution and Dependencies

The final OncoPrint figure can be reproduced by running the visualization script once the files in `data/processed/` are available.

