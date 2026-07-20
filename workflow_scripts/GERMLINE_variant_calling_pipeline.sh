#!/bin/bash

# Usage: ./germline_vc_pipeline.sh <sample_id>
# Example: ./germline_vc_pipeline.sh SRR1653259

set -euo pipefail

#=======================================
# Inputs
#=======================================
SAMPLE=$1

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
fastqc $FASTQ1 $FASTQ2 -o $OUTDIR/qc_reports

#=======================================
# TRIMMING
#=======================================
echo "Trimming the reads..."
trimmomatic PE -threads 4 \
  $FASTQ1 $FASTQ2 \
  $OUTDIR/trimmed_reads/${SAMPLE}_R1_paired.fq.gz $OUTDIR/trimmed_reads/${SAMPLE}_R1_unpaired.fq.gz \
  $OUTDIR/trimmed_reads/${SAMPLE}_R2_paired.fq.gz $OUTDIR/trimmed_reads/${SAMPLE}_R2_unpaired.fq.gz \
  ILLUMINACLIP:$WORKDIR/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 MINLEN:36

#=======================================
# FASTQC on trimmed reads
#=======================================
echo "Running FastQC on trimmed reads........."
for fq in \
  $OUTDIR/trimmed_reads/${SAMPLE}_R1_paired.fq.gz \
  $OUTDIR/trimmed_reads/${SAMPLE}_R2_paired.fq.gz
do
  fastqc -t 8 "$fq" -o $OUTDIR/qc_reports
done

#=======================================
# Aligning with bwa
#=======================================
echo "Aligning reads........"
bwa mem -t 8 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" $REF \
  $OUTDIR/trimmed_reads/${SAMPLE}_R1_paired.fq.gz $OUTDIR/trimmed_reads/${SAMPLE}_R2_paired.fq.gz \
  | samtools sort -o $OUTDIR/bam_files/${SAMPLE}_sorted.bam

samtools index $OUTDIR/bam_files/${SAMPLE}_sorted.bam

#=======================================
# Marking duplicates with gatk
#=======================================
echo "Marking duplicates.............."
gatk MarkDuplicates \
  -I $OUTDIR/bam_files/${SAMPLE}_sorted.bam \
  -O $OUTDIR/bam_files/${SAMPLE}_dedup.bam \
  -M $OUTDIR/logs/${SAMPLE}_metrics.txt

#=======================================
# BQSR
#=======================================
echo "BQSR......................"
gatk BaseRecalibrator \
  -I $OUTDIR/bam_files/${SAMPLE}_dedup.bam -R $REF \
  --known-sites $DBSNP --known-sites $INDELS --known-sites $MILLS \
  -O $OUTDIR/logs/${SAMPLE}_recal.table
  
gatk ApplyBQSR \
  -I $OUTDIR/bam_files/${SAMPLE}_dedup.bam -R $REF \
  --bqsr-recal-file $OUTDIR/logs/${SAMPLE}_recal.table \
  -O $OUTDIR/bam_files/${SAMPLE}_recal.bam

samtools index $OUTDIR/bam_files/${SAMPLE}_recal.bam

#=======================================
# Paths for VCF outputs
#=======================================
RAW_VCF="$OUTDIR/vcf_files/${SAMPLE}_raw.vcf.gz"
GENOTYPED_VCF="$OUTDIR/vcf_files/${SAMPLE}_genotyped.vcf"
FILTERED_SNP_VCF="$OUTDIR/vcf_files/${SAMPLE}_filtered_snps.vcf"
FILTERED_INDEL_VCF="$OUTDIR/vcf_files/${SAMPLE}_filtered_indels.vcf"
MERGED_FILTERED_VCF="$OUTDIR/vcf_files/${SAMPLE}_filtered_merged.vcf.gz"
ANNOTATED_VCF="$OUTDIR/vcf_files/${SAMPLE}_annotated.vcf"

#====================================================
# Calling germline variants with gatk-HaplotypeCaller
#====================================================

echo "Germline variant calling with HaplotypeCaller.........."
RAW_VCF="$OUTDIR/vcf_files/${SAMPLE}_raw.vcf.gz"
gatk --java-options "-Xmx14g" HaplotypeCaller \
  -R $REF \
  -I $OUTDIR/bam_files/${SAMPLE}_recal.bam \
  -O $RAW_VCF \
  -ERC GVCF
  
#====================================================
# GENOTYPING RAW VCF
#====================================================

echo "Genotyping raw VCF..."
gatk GenotypeGVCFs \
  -R $REF \
  -V $RAW_VCF \
  -O $GENOTYPED_VCF
  
#====================================================
# FILTERING SNP
#====================================================
echo "FILTERING SNPs........."
gatk SelectVariants \
  -R $REF \
  -V $GENOTYPED_VCF \
  --select-type-to-include SNP \
  -O "$OUTDIR/vcf_files/${SAMPLE}_raw_snps.vcf"

gatk VariantFiltration \
  -R $REF \
  -V "$OUTDIR/vcf_files/${SAMPLE}_raw_snps.vcf" \
  --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
  --filter-name "SNP_Filter" \
  -O $FILTERED_SNP_VCF

#====================================================
# FILTERING INDELS
#====================================================
echo "FILTERING INDELS........."
gatk SelectVariants \
  -R $REF \
  -V $GENOTYPED_VCF \
  --select-type-to-include INDEL \
  -O "$OUTDIR/vcf_files/${SAMPLE}_raw_indels.vcf"

gatk VariantFiltration \
  -R $REF \
  -V "$OUTDIR/vcf_files/${SAMPLE}_raw_indels.vcf" \
  --filter-expression "QD < 2.0 || FS > 200.0" \
  --filter-name "INDEL_Filter" \
  -O $FILTERED_INDEL_VCF

#====================================================
# MERGING FILTERED SNPS AND INDELS
#====================================================
echo "Merging filtered SNPs and INDELs..."
gatk MergeVcfs \
  -I $FILTERED_SNP_VCF \
  -I $FILTERED_INDEL_VCF \
  -O $MERGED_FILTERED_VCF

# Index merged VCF
if [ ! -f "${MERGED_FILTERED_VCF}.tbi" ]; then
    tabix -p vcf $MERGED_FILTERED_VCF
fi

#====================================================
# ANNOTATING VARIANTS WITH FUNCOTATOR
#====================================================
echo "Annotating filtered variants with Funcotator....."
gatk Funcotator \
  --variant $MERGED_FILTERED_VCF \
  --reference $REF \
  --ref-version hg38 \
  --data-sources-path $FUNCOTATOR \
  --output $ANNOTATED_VCF \
  --output-file-format VCF
  
echo "Germline Variant Calling Pipeline was completed successfully.....!" 
