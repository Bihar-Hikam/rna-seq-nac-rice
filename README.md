# RNA-Seq Analysis of NAC Transcription Factor Genes in Rice Under Drought and Salinity Stress

## Overview
This project presents an end-to-end RNA-Seq analysis to identify differentially expressed NAC transcription factor genes in rice (Oryza sativa L.) under drought and salinity stress conditions. The analysis integrates transcriptomic profiling, variant analysis, promoter characterization, and phylogenetic inference to uncover molecular mechanisms underlying abiotic stress responses.

## Background
Rice (Oryza sativa L.) is a staple crop globally, but its productivity is significantly affected by abiotic stresses such as drought and salinity. NAC (NAM, ATAF1/2, CUC) transcription factors play critical roles in stress response regulation. This study aims to explore NAC gene expression and regulatory mechanisms in stress-tolerant rice varieties.

## Objectives
- Identify differentially expressed genes (DEGs) under drought and salinity stress
- Detect NAC genes responsive to stress conditions
- Analyze SNPs in gene coding and promoter regions
- Predict cis-regulatory elements in promoter regions
- Construct phylogenetic relationships of selected NAC genes

## Data Source
- Public RNA-Seq datasets retrieved from NCBI Sequence Read Archive (SRA)
- Rice varieties: Pokkali (salt-tolerant) and N22 (drought-tolerant)

## Pipeline
1. **Data Acquisition**
   - Download raw sequencing data (FASTQ) from SRA

2. **Quality Control**
   - Tool: FastQC
   - Output: sequencing quality reports

3. **Alignment**
   - Tool: HISAT2
   - Output: SAM/BAM files

4. **Quantification**
   - Tool: HTSeq
   - Output: gene count matrix

5. **Differential Expression Analysis**
   - Tool: R (DESeq2)
   - Output: DEG list (upregulated & downregulated genes)

6. **Functional & Downstream Analysis**
   - SNP analysis (gene & promoter regions)
   - Cis-regulatory element prediction (New PLACE)
   - Protein structure and visualization (PyMOL, AlphaFold)
   - Phylogenetic analysis (MEGA 11, iTOL)

## Tools & Technologies
- Bash (data processing)
- R (DESeq2 for statistical analysis)
- FastQC (quality control)
- HISAT2 / STAR (alignment)
- SAMtools (file processing)
- featureCounts (quantification)
- IGV (visualization)
- PyMOL, AlphaFold (protein analysis)
- MEGA 11, iTOL (phylogenetics)
- New PLACE (cis-element analysis)

## Key Results
- Identified 8,891 DEGs in N22 and 2,081 DEGs in Pokkali
- Four NAC genes consistently expressed across both varieties:
  - LOC_Os01g60020
  - LOC_Os03g60080
  - LOC_Os07g12340
  - LOC_Os11g05614
- LOC_Os03g60080 selected as key candidate due to consistent upregulation and nuclear localization
- SNP analysis revealed:
  - 1 synonymous mutation in coding region
  - 12 SNPs in promoter region potentially affecting gene regulation
- Identified cis-regulatory elements associated with stress response:
  - DPBFCOREDCDC3, GT1GMSCAM4, MYCCONSENSUSAT, RAV1AAT, MYCATERD1
- Phylogenetic analysis showed close relationship with stress-related NAC genes in rice and Arabidopsis

## Key Insights
This study highlights the role of NAC transcription factors, particularly LOC_Os03g60080, in mediating rice responses to drought and salinity stress. The integration of transcriptomic and regulatory analysis provides insights into candidate genes for stress-resilient crop development.

## Notes
- Raw sequencing data is not included due to size limitations.
- Data can be accessed via NCBI SRA (accession IDs provided in scripts).
