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
