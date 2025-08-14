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
