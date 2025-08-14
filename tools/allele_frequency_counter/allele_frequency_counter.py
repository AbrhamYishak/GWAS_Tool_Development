import pandas as pd
import argparse
import sys

parser = argparse.ArgumentParser(description="Allele Frequency Calculator")
parser.add_argument("--genotypes", required=True, help="Path to genotypes.tsv")
parser.add_argument("--snp_annotation", required=False, help="Path to snp_annotation.tsv")
parser.add_argument("--sample_subset", required=False, help="Path to sample subset file (one sample ID per line)")
parser.add_argument("--missing_value", default="NA", help="Missing value code")
parser.add_argument("--output_file", required=True, help="Output TSV file path")
args = parser.parse_args()
if args.snp_annotation in [None, "None"]:
    args.snp_annotation = None
if args.sample_subset in [None, "None"]:
    args.sample_subset = None
try:
    geno = pd.read_csv(args.genotypes, sep="\t", index_col=0, na_values=args.missing_value)
except Exception as e:
    sys.exit(f"Error reading genotypes file: {e}")

snp_annot = None
if args.snp_annotation:
    try:
        snp_annot = pd.read_csv(args.snp_annotation, sep="\t")
        if 'snp_id' not in snp_annot.columns:
            sys.exit("Error: SNP annotation file must contain 'snp_id' column")
    except Exception as e:
        sys.exit(f"Error reading SNP annotation file: {e}")

if args.sample_subset:
    try:
        subset_samples = pd.read_csv(args.sample_subset, header=None)[0].tolist()
        geno = geno.loc[geno.index.isin(subset_samples)]
    except Exception as e:
        sys.exit(f"Error reading sample subset file: {e}")
results = []

for snp in geno.columns:
    g = geno[snp].dropna()
    n_non_missing = g.shape[0]
    if n_non_missing < 10:
        maf = "NA"
    else:
        count_0 = (g == 0).sum()
        count_1 = (g == 1).sum()
        count_2 = (g == 2).sum()
        total_alleles = 2 * n_non_missing
        maf_val = (count_1 + 2*count_2) / total_alleles
        maf = min(maf_val, 1-maf_val)
    
    chrom = snp_annot.loc[snp_annot['snp_id'] == snp, 'chrom'].values[0] if snp_annot is not None else "NA"
    pos = snp_annot.loc[snp_annot['snp_id'] == snp, 'pos'].values[0] if snp_annot is not None else "NA"
    
    results.append({
        "snp_id": snp,
        "chrom": chrom,
        "pos": pos,
        "count_0": count_0,
        "count_1": count_1,
        "count_2": count_2,
        "n_non_missing": n_non_missing,
        "maf": maf
    })
df_out = pd.DataFrame(results)
df_out.to_csv(args.output_file, sep="\t", index=False)

print(f"Allele frequency calculation completed. Output saved to {args.output_file}")
