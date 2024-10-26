#!/bin/bash
# WES Variant Calling Pipeline - GATK Best Practices

# Check if the conda environment exists
ENV_NAME="wes_pipeline_env"
if conda env list | grep -q $ENV_NAME; then
    echo "Activating existing conda environment: $ENV_NAME"
    conda activate $ENV_NAME
else
    echo "Creating new conda environment: $ENV_NAME"
    conda create -n $ENV_NAME -c bioconda gatk4 bwa samtools fastqc picard bcftools vcf-validator -y
    conda activate $ENV_NAME
fi

set -e

function error_exit {
    echo "Error: $1"
    exit 1
}

# check number of available CPU threads
THREADS=$(nproc)
echo "Using $THREADS threads for parallel processes."

HOME_DIR="$HOME"
WES_DIR="$HOME_DIR/wes"
READS_DIR="$HOME_DIR/reads"
ALIGNED_READS_DIR="$HOME_DIR/aligned_reads"
RESULTS_DIR="$HOME_DIR/results"
SUPPORTING_FILES_DIR="$HOME_DIR/supporting_files/hg38"

mkdir -p ${WES_DIR} ${READS_DIR} ${ALIGNED_READS_DIR} ${RESULTS_DIR} ${SUPPORTING_FILES_DIR}

cd $HOME_DIR
echo "pwd: $(pwd)"

REF_GENOME="$HOME_DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

TSV_FILE="$HOME_DIR/igsr_HG00479.tsv"
fastq_1_url=$(awk 'NR==2 {print $1}' $TSV_FILE)
fastq_2_url=$(awk 'NR==3 {print $1}' $TSV_FILE)
sample_name=$(awk 'NR==2 {print $9}' $TSV_FILE)

echo "Sample Name: $sample_name"
echo "FASTQ URLs: $fastq_1_url and $fastq_2_url"

echo "Downloading FASTQ files..."
if [[ ! -f ${READS_DIR}/${sample_name}_1.fastq.gz ]]; then
    echo "Downloading FASTQ file 1..."
    wget -O ${READS_DIR}/${sample_name}_1.fastq.gz ${fastq_1_url} || error_exit "Failed to download FASTQ file 1."
fi

if [[ ! -f ${READS_DIR}/${sample_name}_2.fastq.gz ]]; then
    echo "Downloading FASTQ file 2..."
    wget -O ${READS_DIR}/${sample_name}_2.fastq.gz ${fastq_2_url} || error_exit "Failed to download FASTQ file 2."
fi
echo "Download complete."

if [[ ! -f ${REF_GENOME} ]]; then
    echo "Downloading reference genome..."
    wget -O ${REF_GENOME}.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz || error_exit "Failed to download reference genome."
    gunzip ${REF_GENOME}.gz || error_exit "Failed to unzip reference genome."
fi

# Download BWA index if it does not exist
if [[ ! -f "${REF_GENOME}.bwt" ]]; then
    echo "Downloading BWA index..."
    wget -O ${REF_GENOME}.bwa_index.tar.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz || error_exit "Failed to download BWA index."
    tar -xzvf ${REF_GENOME}.bwa_index.tar.gz -C $(dirname ${REF_GENOME}) || error_exit "Failed to extract BWA index."
    echo "BWA index download complete."
fi

# Generate the FASTA index (.fai) and dictionary (.dict) for GATK
if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "Creating FASTA index..."
    samtools faidx ${REF_GENOME}
fi

if [ ! -f "${REF_GENOME%.fna}.dict" ]; then
    echo "Creating reference dictionary..."
    gatk CreateSequenceDictionary -R ${REF_GENOME} -O ${REF_GENOME%.fna}.dict
fi

# Download known sites for BQSR
KNOWN_SITES="${SUPPORTING_FILES_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf"

if [[ ! -f "${KNOWN_SITES}" ]]; then
    echo "Downloading dbSNP for BQSR..."
    wget -P ${SUPPORTING_FILES_DIR} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf || error_exit "Failed to download dbSNP VCF."
    wget -P ${SUPPORTING_FILES_DIR} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx || error_exit "Failed to download dbSNP VCF."
    echo "dbSNP VCF download complete."
fi


echo "Running FastQC on FASTQ files..."
fastqc ${READS_DIR}/${sample_name}_1.fastq.gz -o ${READS_DIR} || error_exit "FastQC failed for read 1."
fastqc ${READS_DIR}/${sample_name}_2.fastq.gz -o ${READS_DIR} || error_exit "FastQC failed for read 2"
echo "FastQC complete."

echo "Aligning reads to reference genome using BWA-MEM..."
bwa mem -t ${THREADS} -R "@RG\tID:${sample_name}\tPL:ILLUMINA\tSM:${sample_name}" ${REF_GENOME} ${READS_DIR}/${sample_name}_1.fastq.gz ${READS_DIR}/${sample_name}_2.fastq.gz > ${ALIGNED_READS_DIR}/${sample_name}.paired.sam || error_exit "BWA-MEM alignment failed."
echo "Alignment complete."


echo "Marking duplicates..."
gatk MarkDuplicatesSpark -I ${ALIGNED_READS_DIR}/${sample_name}.paired.sam -O ${ALIGNED_READS_DIR}/${sample_name}_sorted_dedup_reads.bam || error_exit "MarkDuplicates failed."
echo "MarkDuplicates complete."


# Remove SAM file to save space
rm ${ALIGNED_READS_DIR}/${sample_name}.paired.sam
echo "SAM file removed."

# Collect alignment metrics
echo "Collecting alignment summary metrics..."
gatk CollectAlignmentSummaryMetrics -I ${ALIGNED_READS_DIR}/${sample_name}_sorted_dedup_reads.bam \
    -R ${REF_GENOME} \
    -O ${RESULTS_DIR}/alignment_summary_metrics.txt || error_exit "Alignment summary metrics collection failed."
echo "Alignment summary metrics collected."

# Base Quality Score Recalibration (BQSR)
echo "Running BaseRecalibrator..."
gatk BaseRecalibrator -I ${ALIGNED_READS_DIR}/${sample_name}_sorted_dedup_reads.bam \
    -R ${REF_GENOME} \
    --known-sites ${KNOWN_SITES} \
     -O recal_data.table || error_exit "BaseRecalibrator failed."

echo "Applying BQSR..."
gatk ApplyBQSR -I ${ALIGNED_READS_DIR}/${sample_name}_sorted_dedup_reads.bam  \
    -R ${REF_GENOME} --bqsr-recal-file recal_data.table \
    -O ${ALIGNED_READS_DIR}/${sample_name}_sorted_dedup_bqsr_reads.bam || error_exit "BQSR failed."
echo "BQSR complete."

# Variant Calling with GATK HaplotypeCaller
echo "Calling variants with GATK HaplotypeCaller..."
gatk HaplotypeCaller -R ${REF_GENOME} \
    -I ${ALIGNED_READS_DIR}/${sample_name}_sorted_dedup_bqsr_reads.bam \
    --dbsnp ${KNOWN_SITES} \
    -O ${RESULTS_DIR}/raw_variants.vcf || error_exit "HaplotypeCaller failed."
echo "Variant calling complete."

# Extract SNPs and Indels
echo "Extracting SNPs and Indels..."
echo "Extracting SNPs..."
gatk SelectVariants -R ${REF_GENOME} \
    -V ${RESULTS_DIR}/raw_variants.vcf  \
    --select-type SNP \
    -O ${RESULTS_DIR}/raw_snps.vcf || error_exit "Failed to extract SNPs."

echo "Extracting INDELs..."
gatk SelectVariants -R ${REF_GENOME} \
    -V ${RESULTS_DIR}/raw_variants.vcf \
    --select-type INDEL \
    -O ${RESULTS_DIR}/raw_indels.vcf || error_exit "Failed to extract Indels."
echo "SNPs and Indels extraction complete."

# Filter SNPs
echo "Filtering SNPs..."
gatk VariantFiltration \
    -R ${REF_GENOME} \
    -V ${RESULTS_DIR}/raw_snps.vcf \
    -O ${RESULTS_DIR}/filtered_snps.vcf \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 60.0" \
    -filter-name "MQ_filter" -filter "MQ < 40.0" \
    -filter-name "SOR_filter" -filter "SOR > 4.0" \
    -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	  -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
    --genotype-filter-expression "DP < 10" \
    --genotype-filter-name "DP_filter" \
    --genotype-filter-expression "GQ < 10" \
    --genotype-filter-name "GQ_filter" || error_exit "Failed to filter SNPs."
echo "SNP filtering complete."


# Filter Indels
echo "Filtering Indels..."
gatk VariantFiltration \
    -R ${REF_GENOME} \
    -V ${RESULTS_DIR}/raw_indels.vcf \
    -O ${RESULTS_DIR}/filtered_indels.vcf \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 200.0" \
    -filter-name "SOR_filter" -filter "SOR > 10.0" \
    --genotype-filter-expression "DP < 10" \
    --genotype-filter-name "DP_filter" \
    --genotype-filter-expression "GQ < 10" \
    --genotype-filter-name "GQ_filter" || error_exit "Failed to filter Indels."
echo "Indel filtering complete."

# Select variants passed filters
echo "Selecting variants passed filters..."
gatk SelectVariants --exclude-filtered -V ${RESULTS_DIR}/filtered_snps.vcf \
    -O ${RESULTS_DIR}/analysis_ready_snps.vcf || error_exit "Failed to select SNPs."
gatk SelectVariants --exclude-filtered -V ${RESULTS_DIR}/filtered_indels.vcf \
    -O ${RESULTS_DIR}/analysis_ready_indels.vcf || error_exit "Failed to select Indels."
echo "Variant selection complete."

# Filter failed SNPs and Indels
echo "Filtering failed variants..."
grep -v -E "DP_filter|GQ_filter" ${RESULTS_DIR}/analysis_ready_snps.vcf > ${RESULTS_DIR}/analysis_ready_snps_filteredGT.vcf  || error_exit "Failed to filter SNPs."
grep -v -E "DP_filter|GQ_filter" ${RESULTS_DIR}/analysis_ready_indels.vcf > ${RESULTS_DIR}/analysis_ready_indels_filteredGT.vcf || error_exit "Failed to filter Indels."
echo "Variant filtering complete."

# Combine filtered SNPs and Indels into a single VCF
echo "Combining filtered SNPs and Indels..."
bcftools sort -Oz ${RESULTS_DIR}/analysis_ready_snps_filteredGT.vcf -o ${RESULTS_DIR}/sorted_analysis_ready_snps_filteredGT.vcf.gz || error_exit "Failed to sort SNPs VCF."
tabix -p vcf "${RESULTS_DIR}/sorted_analysis_ready_snps_filteredGT.vcf.gz" || error_exit "Failed to index SNPs VCF."

bcftools sort -Oz ${RESULTS_DIR}/analysis_ready_indels_filteredGT.vcf -o ${RESULTS_DIR}/sorted_analysis_ready_indels_filteredGT.vcf.gz || error_exit "Failed to sort Indels VCF."
tabix -p vcf "${RESULTS_DIR}/sorted_analysis_ready_indels_filteredGT.vcf.gz" || error_exit "Failed to index Indels VCF."

bcftools concat -a -Oz \
    ${RESULTS_DIR}/sorted_analysis_ready_snps_filteredGT.vcf.gz \
    ${RESULTS_DIR}/sorted_analysis_ready_indels_filteredGT.vcf.gz \
    -o ${RESULTS_DIR}/combined_analysis_ready_${sample_name}.vcf.gz || error_exit "Failed to concatenate VCF files."
echo "Combining filtered SNPs and Indels complete."

tabix -p vcf ${RESULTS_DIR}/combined_analysis_ready_${sample_name}.vcf.gz || error_exit "Failed to index combined VCF."
echo "VCF file ready for upload to wANNOVAR: ${RESULTS_DIR}/combined_analysis_ready_${sample_name}.vcf.gz"

echo "Validating the final VCF file..."
gatk ValidateVariants \
    -V ${RESULTS_DIR}/combined_analysis_ready_${sample_name}.vcf.gz \
    -R ${REF_GENOME} || error_exit "VCF validation failed."
    
## idk why it's not working - need to check on this    
# vcf_validator ${RESULTS_DIR}/combined_analysis_ready_${sample_name}.vcf.gz || error_exit "VCF validation failed."

echo "VCF validation successful."

echo "WES pipeline completed successfully."

conda deactivate
