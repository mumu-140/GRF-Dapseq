#!/bin/bash

# ------------------------------ Step 0: Initialize Parameters and Directory Structure ------------------------------
# Set working directory
WORKDIR=$(pwd)
# Script directory
Script_DIR=${WORKDIR}/scripts
# Sample identifier
name="SampleX"  # Replace with the actual sample prefix

# Raw sequencing data
READ1_REP1="${WORKDIR}/raw/${name}_R1_rep1.fq.gz"
READ2_REP1="${WORKDIR}/raw/${name}_R2_rep1.fq.gz"
READ1_REP2="${WORKDIR}/raw/${name}_R1_rep2.fq.gz"
READ2_REP2="${WORKDIR}/raw/${name}_R2_rep2.fq.gz"
INPUT1="${WORKDIR}/raw/${name}_Input_R1.fq.gz"
INPUT2="${WORKDIR}/raw/${name}_Input_R2.fq.gz"

# Reference genome and annotation files
GENOME_FASTA="${WORKDIR}/reference/genome.fa"
GFF_FILE="${WORKDIR}/reference/genes.gff"

# Script and output directories
OUTPUT_DIR="${WORKDIR}/Dapseq_out/${name}"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit

# Number of threads for parallel processing
THREADS=14
# annotation script
Anno_peak="${Script_DIR}/Anno_peaks.R"
# promoter settings
upst=2000
dost=2000

# ------------------------------ Step 1: Quality Control and Adapter Trimming ------------------------------

echo "[Step 1] Quality control and adapter trimming"

# Perform adapter trimming and filtering using fastp for each replicate and input control
fastp -i $READ1_REP1 -I $READ2_REP1 -o ${name}-1_1_cl.fq.gz -O ${name}-1_2_cl.fq.gz \
      --json ${name}-1.json --html ${name}-1.html --thread $THREADS

fastp -i $READ1_REP2 -I $READ2_REP2 -o ${name}-2_1_cl.fq.gz -O ${name}-2_2_cl.fq.gz \
      --json ${name}-2.json --html ${name}-2.html --thread $THREADS

fastp -i $INPUT1 -I $INPUT2 -o ${name}_Input_1_cl.fq.gz -O ${name}_Input_2_cl.fq.gz \
      --json ${name}_Input.json --html ${name}_Input.html --thread $THREADS


# ------------------------------ Step 2: Read Alignment to the Reference Genome ------------------------------

echo "[Step 2] Read alignment using BWA"

# Build BWA index for the reference genome
bwa index -p "genome" "$GENOME_FASTA"

# Align cleaned reads to the reference genome using BWA-MEM
bwa mem -t $THREADS -M "genome" ${name}-1_1_cl.fq.gz ${name}-1_2_cl.fq.gz > ${name}-1.sam
bwa mem -t $THREADS -M "genome" ${name}-2_1_cl.fq.gz ${name}-2_2_cl.fq.gz > ${name}-2.sam

# Align input control reads
bwa mem -t $THREADS -M "genome" ${name}_Input_1_cl.fq.gz ${name}_Input_2_cl.fq.gz > ${name}_Input.sam


# ------------------------------ Step 3: BAM Conversion, Sorting, and Duplicate Removal ------------------------------

echo "[Step 3] SAM to BAM conversion, sorting, and duplicate removal"

# Convert SAM to sorted BAM and remove PCR duplicates using Picard tools
for sample in 1 2 Input; do
    picard SortSam SO=coordinate INPUT=${name}-${sample}.sam OUTPUT=${name}-${sample}.sorted.bam \
        VALIDATION_STRINGENCY=LENIENT  CREATE_INDEX=true

    picard MarkDuplicates I=${name}-${sample}.sorted.bam O=${name}-${sample}.dedup.bam M=${name}-${sample}.dedup.txt \
        VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true  CREATE_INDEX=true
    samtools index ${name}-${sample}.dedup.bam
done


# ------------------------------ Step 4: Peak Calling using MACS2 ------------------------------

echo "[Step 4] Peak calling with MACS2"

# Estimate effective genome size using seqkit
GENOME_SIZE=$(seqkit stats "$GENOME_FASTA" | awk 'NR==2 {gsub(/,/, "", $5); print $5}')

# Perform peak calling using MACS2 with paired-end BAM files and input control
macs2 callpeak -f BAMPE -t ${name}-1.dedup.bam -c ${name}-Input.dedup.bam -g $GENOME_SIZE \
    -n ${name}-1 --outdir ./ -B -q 0.05

macs2 callpeak -f BAMPE -t ${name}-2.dedup.bam -c ${name}-Input.dedup.bam -g $GENOME_SIZE \
    -n ${name}-2 --outdir ./ -B -q 0.05

echo "
===================== narrowPeak =====================
Output peak file format: narrowPeak (BED6+4)
    Column 1: Chromosome
    Column 2: Peak start
    Column 3: Peak end
    Column 4: Peak ID
    Column 5: Score (0-1000)
    Column 6: Strand (+/- or '.')
    Column 7: Signal value
    Column 8: -log10(p-value)
    Column 9: -log10(q-value)
    Column 10: Relative summit position
==============================================================
"
# ------------------------------ Step 5: Reproducibility Assessment using IDR ------------------------------

echo "[Step 5] IDR analysis for reproducibility between replicates"

# Sort narrowPeak files by signal p-value (column 8) in descending order
sort -k8,8nr ./${name}-1_peaks.narrowPeak > ${name}-1_sorted.narrowPeak
sort -k8,8nr ./${name}-2_peaks.narrowPeak > ${name}-2_sorted.narrowPeak

# Perform IDR analysis to assess reproducibility of peaks between replicates
idr --samples ${name}-1_sorted.narrowPeak ${name}-2_sorted.narrowPeak \
    --input-file-type narrowPeak --rank p.value \
    --output-file ${name}_idr --plot --log-output-file ${name}_idr.log

echo "IDR analysis completed successfully"

# ------------------------------ Step 6: Peak Annotation ------------------------------
Rscript ${Anno_peak} $GFF_FILE ${name}_idr $OUTPUT_DIR $upst $dost

echo "Annotation of peaks completed successfully"