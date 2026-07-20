#!/bin/bash

# Usage: ./somatic_vc_pipeline.sh <tumor_id> <normal_id>
# Example: ./somatic_vc_pipeline.sh SRR1653259 SRR1656608

set -euo pipefail

#=======================================
# Inputs
#=======================================
TUMOR=SRR1653259
NORMAL=SRR1656608

#=======================================
# Defining paths
#=======================================
WORKDIR="/home/oindree-bal/WORKDIR/Samples/Tumor"
REF="/home/oindree-bal/WORKDIR/Samples/Tumor/REF/Homo_sapiens_assembly38.fasta"
DBSNP="/home/oindree-bal/WORKDIR/DBSNP/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
INDELS="/home/oindree-bal/WORKDIR/INDELS/Homo_sapiens_assembly38.known_indels.vcf.gz"
MILLS="/home/oindree-bal/WORKDIR/MILLS/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
FUNCOTATOR="/home/oindree-bal/WORKDIR/funcotator_dataSources.v1.8.hg38.20230908g"
PANEL="/home/oindree-bal/WORKDIR/targets.bed"

cd "$WORKDIR"

#=======================================
# FASTQ naming convention
#=======================================
FASTQ_TUMOR1="${TUMOR}_1.fastq.gz"
FASTQ_TUMOR2="${TUMOR}_2.fastq.gz"
FASTQ_NORMAL1="${NORMAL}_1.fastq.gz"
FASTQ_NORMAL2="${NORMAL}_2.fastq.gz"

#=======================================
# Output structure
#=======================================
OUTPUT="$WORKDIR/outputs"
mkdir -p $OUTPUT/{qc_reports,trimmed_reads,bam_files,vcf_files,logs}

#=======================================================================================================
# WORKFLOW
#=======================================================================================================

echo "Checking required WES tools..."

# List required tools
REQUIRED_TOOLS=("fastp" "trimmomatic" "bwa" "samtools" "gatk")

for tool in "${REQUIRED_TOOLS[@]}"
do
    if command -v "$tool" &> /dev/null
    then
        echo "[OK] $tool found"
    else
        echo "[ERROR] $tool not found. Please install it or add it to PATH."
        exit 1
    fi
done

echo "All required tools are available."

#=======================================
# FASTQC
#=======================================
echo "FastQC is running.........."
fastqc $FASTQ_TUMOR1 $FASTQ_TUMOR2 $FASTQ_NORMAL1 $FASTQ_NORMAL2 -o $OUTDIR/qc_reports

#=======================================
# TRIMMING
#=======================================
echo "Trimming the reads..."
trimmomatic PE -threads 4 \
  $FASTQ_TUMOR1 $FASTQ_TUMOR2 \
  $OUTPUT/trimmed_reads/${TUMOR}_R1_paired.fq.gz $OUTPUT/trimmed_reads/${TUMOR}_R1_unpaired.fq.gz \
  $OUTPUT/trimmed_reads/${TUMOR}_R2_paired.fq.gz $OUTPUT/trimmed_reads/${TUMOR}_R2_unpaired.fq.gz \
  ILLUMINACLIP:$WORKDIR/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36

trimmomatic PE -threads 4 \
  $FASTQ_NORMAL1 $FASTQ_NORMAL2 \
  $OUTPUT/trimmed_reads/${NORMAL}_R1_paired.fq.gz $OUTPUT/trimmed_reads/${NORMAL}_R1_unpaired.fq.gz \
  $OUTPUT/trimmed_reads/${NORMAL}_R2_paired.fq.gz $OUTPUT/trimmed_reads/${NORMAL}_R2_unpaired.fq.gz \
  ILLUMINACLIP:$WORKDIR/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
  
#=======================================
# FASTQC on trimmed reads
#=======================================
echo "Running FastQC on trimmed reads........."
for fq in \
  $OUTPUT/trimmed_reads/${TUMOR}_R1_paired.fq.gz \
  $OUTPUT/trimmed_reads/${TUMOR}_R2_paired.fq.gz \
  $OUTPUT/trimmed_reads/${NORMAL}_R1_paired.fq.gz \
  $OUTPUT/trimmed_reads/${NORMAL}_R2_paired.fq.gz
do
  fastqc -t 8 "$fq" -o $OUTPUT/qc_reports
done

#=======================================
# Aligning with bwa
#=======================================
echo "Aligning reads........"
bwa mem -t 8 -R "@RG\tID:${TUMOR}\tSM:${TUMOR}\tPL:ILLUMINA" $REF \
  $OUTPUT/trimmed_reads/${TUMOR}_R1_paired.fq.gz $OUTDIR/trimmed_reads/${TUMOR}_R2_paired.fq.gz \
  | samtools sort -o $OUTPUT/bam_files/${TUMOR}_sorted.bam

bwa mem -t 8 -R "@RG\tID:${NORMAL}\tSM:${NORMAL}\tPL:ILLUMINA" $REF \
  $OUTPUT/trimmed_reads/${NORMAL}_R1_paired.fq.gz $OUTDIR/trimmed_reads/${NORMAL}_R2_paired.fq.gz \
  | samtools sort -o $OUTDIR/bam_files/${NORMAL}_sorted.bam
  
#=======================================
# Indexing with samtools
#=======================================
echo "Indexing output bams........."
samtools index $OUTPUT/bam_files/${TUMOR}_sorted.bam
samtools index $OUTPUT/bam_files/${NORMAL}_sorted.bam

#=======================================
# Marking duplicates with gatk
#=======================================
echo "Marking duplicates.............."
gatk MarkDuplicates \
  -I $OUTPUT/bam_files/${TUMOR}_sorted.bam \
  -O $OUTPUT/bam_files/${TUMOR}_dedup.bam \
  -M $OUTPUT/logs/${TUMOR}_metrics.txt

gatk MarkDuplicates\
  -I $OUTPUT/bam_files/${NORMAL}_sorted.bam \
  -O $OUTPUT/bam_files/${NORMAL}_dedup.bam \
  -M $OUTPUT/logs/${NORMAL}_metrics.txt
  
#=======================================
# BQSR
#=======================================

echo "BQSR......................"
gatk BaseRecalibrator \
  -I $OUTPUT/bam_files/${TUMOR}_dedup.bam -R $REF \
  --known-sites $DBSNP --known-sites $INDELS --known-sites $MILLS \
  -O $OUTPUT/logs/${TUMOR}_recal.table
  
gatk ApplyBQSR \
  -I $OUTPUT/bam_files/${TUMOR}_dedup.bam -R $REF \
  --bqsr-recal-file $OUTPUT/logs/${TUMOR}_recal.table \
  -O $OUTPUT/bam_files/${TUMOR}_recal.bam
samtools index $OUTPUT/bam_files/${TUMOR}_recal.bam

gatk BaseRecalibrator \
  -I $OUTPUT/bam_files/${NORMAL}_dedup.bam -R $REF \
  --known-sites $DBSNP --known-sites $INDELS --known-sites $MILLS \
  -O $OUTPUT/logs/${NORMAL}_recal.table
  
gatk ApplyBQSR \
  -I $OUTPUT/bam_files/${NORMAL}_dedup.bam -R $REF \
  --bqsr-recal-file $OUTDIR/logs/${NORMAL}_recal.table \
  -O $OUTPUT/bam_files/${NORMAL}_recal.bam
samtools index $OUTDIR/bam_files/${NORMAL}_recal.bam

#=======================================
# Paths for VCF outputs
#=======================================

RAW_VCF="$OUTPUT/vcf_files/${TUMOR}_vs_${NORMAL}_somatic_raw.vcf.gz"
FILTERED_VCF="$OUTPUT/vcf_files/${TUMOR}_vs_${NORMAL}_somatic_filtered.vcf.gz"
ANNOTATED_VCF="$OUTPUT/vcf_files/${TUMOR}_vs_${NORMAL}_somatic_annotated.vcf.gz"


#=============================================
# Calling somatic variants with gatk-Mutect2
#=============================================

echo "Somatic variant calling .........."
gatk --java-options "-Xmx8g" Mutect2 \
  -R $REF \
  -I $OUTPUT/bam_files/${TUMOR}_recal.bam -tumor $TUMOR \
  -I $OUTPUT/bam_files/${NORMAL}_recal.bam -normal $NORMAL \
  -O $OUTPUT/vcf_files/${TUMOR}_vs_${NORMAL}_somatic_raw.vcf.gz
  
#========================================================
# Filtering somatic variants with gatk-FilterMutectCalls
#========================================================

echo "Filtering Mutect2 calls..."
gatk FilterMutectCalls \
  -R $REF \
  -V $RAW_VCF \
  -O $FILTERED_VCF

#========================================================
# Indexing the variants
#========================================================

if [ ! -f "${FILTERED_VCF}.tbi" ]; then
    tabix -p vcf $FILTERED_VCF
fi

#========================================================
# ANNOTATING with funcotator
#========================================================

echo "Annotating with Funcotator......................"
gatk Funcotator \
  --variant $FILTERED_VCF \
  --reference $REF \
  --ref-version hg38 \
  --data-sources-path $FUNCOTATOR \
  --output $ANNOTATED_VCF \
  --output-file-format VCF
  
echo "Somatic Variant Calling Pipeline was completed successfully.....!" 
