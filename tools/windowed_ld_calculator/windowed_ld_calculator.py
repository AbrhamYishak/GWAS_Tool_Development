import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import sys
parser = argparse.ArgumentParser(description="Windowed LD Calculator")
parser.add_argument("--genotypes", required=True, help="Genotype matrix TSV")
parser.add_argument("--snp_annotation", required=True, help="SNP annotation TSV")
parser.add_argument("--focal_snp", required=True, help="Focal SNP ID or Chrom:Position")
parser.add_argument("--window_kb", type=int, default=250, help="Window size in kb")
parser.add_argument("--min_maf", type=float, default=0.01, help="Minimum MAF")
parser.add_argument("--missing_code", default="NA", help="Missing value code")
parser.add_argument("--ld_matrix_out", required=True, help="Output LD matrix TSV")
parser.add_argument("--heatmap_out", required=True, help="Output heatmap PNG")
args = parser.parse_args()
geno = pd.read_csv(args.genotypes, sep="\t", index_col=0, na_values=args.missing_code)
annot = pd.read_csv(args.snp_annotation, sep="\t")
if ":" in args.focal_snp:
    chrom, pos = args.focal_snp.split(":")
    pos = int(pos)
    focal_snps = annot[(annot["chrom"].astype(str) == str(chrom)) &
                       (annot["pos"] >= pos - args.window_kb*1000) &
                       (annot["pos"] <= pos + args.window_kb*1000)]
else:
    focal_snps = annot[annot["snp_id"] == args.focal_snp]

if focal_snps.empty:
    sys.exit("No SNPs found in the specified window around focal SNP.")
window_snps = annot[(annot["chrom"] == focal_snps["chrom"].values[0]) &
                    (annot["pos"] >= focal_snps["pos"].values[0] - args.window_kb*1000) &
                    (annot["pos"] <= focal_snps["pos"].values[0] + args.window_kb*1000)]

snps_to_use = window_snps["snp_id"].tolist()
geno_window = geno[snps_to_use]
maf = geno_window.apply(lambda x: x.dropna().sum() / (2*x.dropna().shape[0]))
geno_window = geno_window.loc[:, maf >= args.min_maf]
def compute_r2(x, y):
    valid = x.notna() & y.notna()
    x = x[valid]
    y = y[valid]
    if len(x) < 3:
        return np.nan
    cov = ((x - x.mean()) * (y - y.mean())).mean()
    r2 = cov**2 / (x.var() * y.var()) if x.var() > 0 and y.var() > 0 else np.nan
    return r2
ld_matrix = pd.DataFrame(index=geno_window.columns, columns=geno_window.columns, dtype=float)
for snp1 in geno_window.columns:
    for snp2 in geno_window.columns:
        if snp1 <= snp2:
            r2 = compute_r2(geno_window[snp1], geno_window[snp2])
            ld_matrix.loc[snp1, snp2] = r2
            ld_matrix.loc[snp2, snp1] = r2
ld_matrix.to_csv(args.ld_matrix_out, sep="\t")
plt.figure(figsize=(10, 8))
sns.heatmap(ld_matrix.astype(float), cmap="Reds", vmin=0, vmax=1)
plt.title(f"LD Heatmap around {args.focal_snp}")
plt.tight_layout()
if not args.heatmap_out.lower().endswith(".png"):
    args.heatmap_out += ".png"
plt.savefig(args.heatmap_out, dpi=300, format='png')
plt.close()
print(f"LD matrix saved to {args.ld_matrix_out}")
print(f"LD heatmap saved to {args.heatmap_out}")
