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
