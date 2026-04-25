#!/bin/bash

# RNA-Seq Pipeline - Oryza sativa

# ===== CONFIG =====
BASE_DIR="$HOME/transkriptomik"
INPUT_DIR="$BASE_DIR/input"
QC_DIR="$BASE_DIR/fastqc_output"
TRIM_DIR="$BASE_DIR/trimmed"
MAP_DIR="$BASE_DIR/mapped"
COUNT_DIR="$BASE_DIR/counts"

GENOME_INDEX="$BASE_DIR/genome/Oryza_sativa_index"
GTF_FILE="$BASE_DIR/genome/Oryza_sativa.IRGSP1.0.61.gtf"
ADAPTER="$HOME/miniconda3/envs/cobaskripsi/share/trimmomatic-0.39-2/adapters/TruSeq3-SE.fa"

THREADS=4

# ===== CREATE DIRECTORIES =====
mkdir -p $QC_DIR $TRIM_DIR $MAP_DIR $COUNT_DIR

# ===== SAMPLE LIST =====
SAMPLES=(
"N22_R1_Drought"
"N22_R2_Drought"
"N22_R3_Drought"
"N22_R1_Kontrol"
"N22_R2_Kontrol"
"N22_R3_Kontrol"
"Pok_R1_Salt"
"Pok_R2_Salt"
"Pok_R3_Salt"
"Pok_R1_Kontrol"
"Pok_R2_Kontrol"
"Pok_R3_Kontrol"
)

# ================================
# PIPELINE START
# ================================

for SAMPLE in "${SAMPLES[@]}"
do
    echo "======================================"
    echo "Processing: $SAMPLE"
    echo "======================================"

    INPUT_FASTQ="$INPUT_DIR/${SAMPLE}.fastqsanger.gz"
    TRIMMED_FASTQ="$TRIM_DIR/${SAMPLE}_trimmed.fastq"
    BAM_FILE="$MAP_DIR/${SAMPLE}_aligned_sorted.bam"
    COUNT_FILE="$COUNT_DIR/${SAMPLE}_counts.txt"

    # ===== FASTQC =====
    echo "[1/4] Running FastQC..."
    fastqc $INPUT_FASTQ -o $QC_DIR

    # ===== TRIMMING =====
    echo "[2/4] Trimming..."
    trimmomatic SE -threads $THREADS -phred33 \
        $INPUT_FASTQ \
        $TRIMMED_FASTQ \
        ILLUMINACLIP:$ADAPTER:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:35

    # ===== ALIGNMENT + SORT =====
    echo "[3/4] Aligning with HISAT2..."
    hisat2 -p $THREADS -x $GENOME_INDEX -U $TRIMMED_FASTQ | \
    samtools view -bS - | \
    samtools sort -o $BAM_FILE

    # ===== INDEX BAM =====
    samtools index $BAM_FILE

    # ===== COUNTING =====
    echo "[4/4] Counting reads..."
    htseq-count -f bam -r pos -s no -t exon -i gene_id \
        $BAM_FILE \
        $GTF_FILE \
        > $COUNT_FILE

    echo "Finished: $SAMPLE"
done

echo "ALL SAMPLES PROCESSED SUCCESSFULLY"
