import pandas as pd
import argparse
import sys
import statsmodels.api as sm
parser = argparse.ArgumentParser(description="Logistic Regression Association Test")
parser.add_argument("--genotypes", required=True, help="Path to genotypes.tsv")
parser.add_argument("--phenotypes_covariates", required=True, help="Path to phenotypes/covariates TSV")
parser.add_argument("--snp_annotation", required=False, help="Path to snp_annotation.tsv")
parser.add_argument("--covariates", required=True, help="Comma-separated covariates to include")
parser.add_argument("--min_maf", type=float, default=0.01, help="Minimum MAF to include SNP")
parser.add_argument("--missing_value", default="NA", help="Missing value code")
parser.add_argument("--output_file", required=True, help="Output TSV file path")
args = parser.parse_args()
if args.snp_annotation in [None, "None"]:
    args.snp_annotation = None
try:
    geno = pd.read_csv(args.genotypes, sep="\t", index_col=0, na_values=args.missing_value)
except Exception as e:
    sys.exit(f"Error reading genotypes file: {e}")

try:
    pheno = pd.read_csv(args.phenotypes_covariates, sep="\t", na_values=args.missing_value)
    if 'phenotype' not in pheno.columns:
        sys.exit("Error: phenotype column must exist in phenotypes_covariates.tsv")
except Exception as e:
    sys.exit(f"Error reading phenotypes/covariates file: {e}")
covariates = [c.strip() for c in args.covariates.split(",") if c.strip() != ""]
missing_covs = [c for c in covariates if c not in pheno.columns]
if missing_covs:
    sys.exit(f"Error: The following covariates are missing in phenotype file: {missing_covs}")
snp_annot = None
if args.snp_annotation:
    try:
        snp_annot = pd.read_csv(args.snp_annotation, sep="\t")
        if 'snp_id' not in snp_annot.columns:
            sys.exit("Error: SNP annotation file must contain 'snp_id' column")
    except Exception as e:
        sys.exit(f"Error reading SNP annotation file: {e}")
common_samples = geno.index.intersection(pheno['sample_id'])
geno = geno.loc[common_samples]
pheno = pheno.loc[pheno['sample_id'].isin(common_samples)].set_index('sample_id')
results = []

for snp in geno.columns:
    g = geno[snp]
    g_non_missing = g.dropna()
    n_non_missing = g_non_missing.shape[0]
    if n_non_missing < 10:
        continue

    total_alleles = 2 * n_non_missing
    allele_sum = g_non_missing.sum()
    maf = min(allele_sum / total_alleles, 1 - allele_sum / total_alleles)
    if maf < args.min_maf:
        continue
    df = pheno.join(g.rename("genotype"))
    df = df.dropna(subset=['phenotype', 'genotype'] + covariates)

    X = df[['genotype'] + covariates]
    X = sm.add_constant(X)
    y = df['phenotype']

    try:
        model = sm.Logit(y, X)
        res = model.fit(disp=False)
        beta = res.params['genotype']
        se = res.bse['genotype']
        z = beta / se
        pval = res.pvalues['genotype']
    except Exception:
        beta = se = z = pval = "NA"

    chrom = snp_annot.loc[snp_annot['snp_id'] == snp, 'chrom'].values[0] if snp_annot is not None else "NA"
    pos = snp_annot.loc[snp_annot['snp_id'] == snp, 'pos'].values[0] if snp_annot is not None else "NA"

    results.append({
        "snp_id": snp,
        "chrom": chrom,
        "pos": pos,
        "beta": beta,
        "SE": se,
        "z": z,
        "p_value": pval,
        "MAF": maf
    })
df_out = pd.DataFrame(results)
df_out.to_csv(args.output_file, sep="\t", index=False)
print(f"Association testing completed. Output saved to {args.output_file}")
