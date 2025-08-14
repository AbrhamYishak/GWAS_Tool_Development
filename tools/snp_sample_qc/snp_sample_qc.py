import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description="SNP & Sample QC Tool")
parser.add_argument("--genotypes", required=True, help="Path to genotypes.tsv")
parser.add_argument("--snp_annotation", required=True, help="Path to snp_annotation.tsv")
parser.add_argument("--phenotypes", required=False, help="Path to phenotypes_covariates.tsv")
parser.add_argument("--maf_threshold", type=float, default=0.01, help="MAF threshold")
parser.add_argument("--snp_missingness_threshold", type=float, default=0.05, help="SNP missingness threshold")
parser.add_argument("--sample_missingness_threshold", type=float, default=0.1, help="Sample missingness threshold")
parser.add_argument("--missing_value", default="NA", help="Missing value code")
parser.add_argument("--output_dir", default="qc_output", help="Output directory")
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

geno = pd.read_csv(args.genotypes, sep="\t", index_col=0, na_values=args.missing_value)
snp_annot = pd.read_csv(args.snp_annotation, sep="\t")
phenos = None
if args.phenotypes:
    phenos = pd.read_csv(args.phenotypes, sep="\t")

report = []

snp_missing = geno.isna().mean(axis=0)
filtered_snps = snp_missing[snp_missing <= args.snp_missingness_threshold].index
report.append(f"SNPs before missingness filter: {geno.shape[1]}, after: {len(filtered_snps)}")
geno = geno[filtered_snps]
snp_annot = snp_annot[snp_annot['snp_id'].isin(filtered_snps)]

sample_missing = geno.isna().mean(axis=1)
filtered_samples = sample_missing[sample_missing <= args.sample_missingness_threshold].index
report.append(f"Samples before missingness filter: {geno.shape[0]}, after: {len(filtered_samples)}")
geno = geno.loc[filtered_samples]
if phenos is not None:
    phenos = phenos[phenos['sample_id'].isin(filtered_samples)]

allele_sums = geno.sum(axis=0)
allele_counts = 2 * geno.notna().sum(axis=0)  # fixed line
maf = allele_sums / allele_counts
maf = maf.apply(lambda x: min(x, 1-x))
filtered_maf_snps = maf[maf >= args.maf_threshold].index
report.append(f"SNPs before MAF filter: {geno.shape[1]}, after: {len(filtered_maf_snps)}")
geno = geno[filtered_maf_snps]
snp_annot = snp_annot[snp_annot['snp_id'].isin(filtered_maf_snps)]

geno.to_csv(os.path.join(args.output_dir, "filtered_genotypes.tsv"), sep="\t")
snp_annot.to_csv(os.path.join(args.output_dir, "filtered_snp_annotation.tsv"), sep="\t", index=False)
out_pheno = os.path.join(args.output_dir, "filtered_phenotypes.tsv")
if phenos is not None:
    phenos.to_csv(out_pheno, sep="\t", index=False)
else:
    pd.DataFrame(columns=["sample_id", "phenotype", "age", "sex"]).to_csv(out_pheno, sep="\t", index=False)

with open(os.path.join(args.output_dir, "qc_report.txt"), "w") as f:
    f.write("\n".join(report))

print("QC completed. Outputs saved in:", args.output_dir)
