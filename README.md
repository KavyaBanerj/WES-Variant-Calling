# WES-Variant-Calling

## Overview

The **WES-Variant-Calling** workflow is designed to process human Whole Exome Sequencing (WES) data following [GATK4 best practices for germline variant calling](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels). This pipeline encompasses data downloading, quality control, alignment, duplicate marking, base quality score recalibration (BQSR), variant calling, filtering, and VCF validation.

## Prerequisites

Before running the WES-Variant-Calling pipeline, ensure you have the following:

1. **Reference Files :**
   - **Human Reference Genome (downloaded during execution):** GRCh38 (`GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`)
   - **dbSNP VCF (downloaded during execution):** `Homo_sapiens_assembly38.dbsnp138.vcf`
   - **Sample Metadata TSV:** Contains URLs for FASTQ files and sample information

2. **Software Dependencies:**
   - **Anaconda:** For environment management
   - **Bioinformatics Tools (downloaded during execution):** Installed via Conda (e.g., GATK4, BWA-MEM, Samtools, FastQC, Picard, BCFtools, VCF-validator)

3. **Compute Resources:**
   - Sufficient memory and CPU cores as required by the pipeline steps

## Installation
### 1. Clone the Repository

```bash
git clone https://github.com/KavyaBanerj/WES-Variant-Calling.git
cd WES-Variant-Calling
```

## Usage

### 1. Prepare Sample Metadata

Create a TSV file (e.g., `igsr_HG00479.tsv`) containing the URLs for your FASTQ files and sample information.


### 2. Configure the Pipeline Script

Ensure that the `wes_pipeline.sh` script has the correct paths and configurations based on your directory structure and sample metadata.

### 3. Run the Pipeline

Execute the pipeline script:

```bash
bash ./wes_pipeline.sh
```

**Note:** Ensure that the script has execute permissions. If not, set them using:

```bash
chmod +x ./wes_pipeline.sh
```

## Pipeline Steps

1. **Data Downloading:**
   - Downloads FASTQ files based on URLs provided in the metadata TSV.
   - Downloads and prepares the reference genome and dbSNP VCF.

2. **Quality Control (QC):**
   - Runs FastQC on raw FASTQ files to assess quality metrics.

3. **Alignment:**
   - Aligns reads to the reference genome using BWA-MEM.
   - Includes read group (`@RG`) information for downstream analyses.

4. **Duplicate Marking:**
   - Marks PCR duplicates using GATK's `MarkDuplicatesSpark`.

5. **Base Quality Score Recalibration (BQSR):**
   - Performs BQSR to correct systematic errors in base quality scores.

6. **Variant Calling:**
   - Calls variants using GATK's `HaplotypeCaller`.
   - Generates raw VCF files with variant annotations.

7. **Variant Extraction:**
   - Extracts SNPs and Indels into separate VCF files using `SelectVariants`.

8. **Variant Filtering:**
   - Applies filters based on quality metrics (`QD`, `FS`, `MQ`, `SOR`, `MQRankSum`, `ReadPosRankSum`) using `VariantFiltration`.
   - Applies genotype-level filters (`DP`, `GQ`).

9. **VCF Processing:**
   - Further filters variants using `grep` to remove genotype filter flags.
   - Sorts, compresses, indexes, and concatenates SNP and Indel VCFs using `bcftools`.

10. **VCF Validation:**
    - Validates the final VCF file using GATK's `ValidateVariants`.

## Outputs

The pipeline generates several output files within the `results/` directory:

- **Raw Variants:**
  - `raw_variants.vcf`: Unfiltered variant calls.

- **Filtered Variants:**
  - `filtered_snps.vcf`: SNPs after applying variant filters.
  - `filtered_indels.vcf`: Indels after applying variant filters.

- **Analysis-Ready Variants:**
  - `analysis_ready_snps.vcf`: SNPs passing all filters.
  - `analysis_ready_indels.vcf`: Indels passing all filters.

- **Final Combined VCF:**
  - `combined_analysis_ready_{sample_name}.vcf.gz`: Concatenated and compressed VCF file ready for downstream tools like wANNOVAR.

## Acknowledgements

- Pipeline inspired by [WES-analysis](https://github.com/mariamnawaz1/WES-analysis) and [variant_calling](https://github.com/kpatel427/YouTubeTutorials/blob/main/variant_calling.sh)
