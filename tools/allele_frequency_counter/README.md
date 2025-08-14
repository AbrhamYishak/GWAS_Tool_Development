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
