# GWAS Tool Development – Assignment

## 1. Requirements

1. Python 3.8+
2. pandas
3. numpy
4. scipy
5. matplotlib
6. statsmodels (for logistic regression)
7. seaborn

---

## 2. Data Preparation

Synthetic GWAS data can be generated locally using the following command:

```bash
python generate_gwas_data.py
```

Generated files:

* `genotypes.tsv` → samples × SNPs matrix (0,1,2 genotype calls)
* `phenotypes_covariates.tsv` → phenotype + covariates
* `snp_annotation.tsv` → SNP info (chrom, pos, maf)

---

## 3. Tool Specifications

### Tool 1: SNP & Sample QC

**Input:**

* genotypes.tsv
* snp\_annotation.tsv
* phenotypes\_covariates.tsv (optional)
* MAF threshold (default = 0.01)
* SNP missingness threshold (default = 0.05)
* Sample missingness threshold (default = 0.1)
* Missing value code (default = NA)

**Output:**

* Filtered genotypes
* Filtered SNP annotation
* Filtered sample list
* QC report (txt with summary stats)

**Logic:**

* Remove SNPs with missingness > threshold
* Remove SNPs with MAF < threshold
* Remove samples with missingness > threshold
* Output counts before/after filtering

---

### Tool 2: Allele Frequency Calculator

**Input:**

* genotypes.tsv
* snp\_annotation.tsv (optional)
* Sample subset file (optional)
* Missing value code (default = NA)

**Output:**

* TSV with SNP ID, chrom, pos, counts of each genotype, n\_non\_missing, and MAF

**Logic:**

* Compute allele counts (0, 1, 2)
* Calculate MAF per SNP
* If < 10 non-missing samples, set MAF = NA

---

### Tool 3: Association Test (Logistic regression)

**Input:**

* genotypes.tsv
* phenotypes\_covariates.tsv
* snp\_annotation.tsv (optional)
* Covariates (comma-separated)
* Min MAF (default = 0.01)
* Missing value code

**Output:**

* Association results: SNP ID, chrom, pos, beta, SE, z, p-value, MAF

**Logic:**

* For each SNP, fit logistic regression: phenotype \~ genotype + covariates
* Skip SNPs with MAF < threshold or insufficient samples

---

### Tool 4: Manhattan Plot Generator

**Input:**

* Association results (from Tool 3)
* SNP annotation (optional)
* P-value threshold (default = 5e-8)
* Top N hits (default = 20)

**Output:**

* Manhattan plot PNG
* Top hits TSV

**Logic:**

* Plot -log10(p) by chrom/pos
* Highlight hits below threshold
* Output top N SNPs by p-value

---

### Tool 5: Windowed LD Calculator

**Input:**

* genotypes.tsv
* snp\_annotation.tsv
* Focal SNP ID or chrom + center pos
* Window size in kb (default = 250)
* Min MAF (default = 0.01)
* Missing value code

**Output:**

* LD matrix TSV
* Heatmap PNG

**Logic:**

* Select SNPs within window
* Compute r² for each pair
* Plot heatmap
