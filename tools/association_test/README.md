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
