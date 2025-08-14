import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import sys
parser = argparse.ArgumentParser(description="Manhattan Plot Generator")
parser.add_argument("--assoc_results", required=True, help="Association results TSV")
parser.add_argument("--annotation", required=False, help="SNP annotation TSV (optional)")
parser.add_argument("--pvalue_threshold", type=float, default=5e-8, help="P-value threshold")
parser.add_argument("--top_n", type=int, default=20, help="Number of top hits to output")
parser.add_argument("--plot_out", required=True, help="Output PNG file path")
parser.add_argument("--hits_out", required=True, help="Output top hits TSV file path")
args = parser.parse_args()
df = pd.read_csv(args.assoc_results, sep="\t")
pval_col = next((c for c in df.columns if c.lower() in ["p_value", "pval", "p"]), None)
if pval_col is None:
    sys.exit("No p-value column found in association results.")
if args.annotation:
    annot = pd.read_csv(args.annotation, sep="\t")
    chrom_col = next((c for c in annot.columns if c.lower() in ["chrom", "chr", "chromosome"]), None)
    pos_col = next((c for c in annot.columns if c.lower() in ["pos", "bp", "position"]), None)
    snp_col = next((c for c in df.columns if c.lower() in annot.columns or c.lower() in ["snp", "marker", "id"]), None)
    if chrom_col and pos_col and snp_col:
        df = df.merge(annot[[snp_col, chrom_col, pos_col]], left_on=snp_col, right_on=snp_col, how="left")
        df.rename(columns={chrom_col: "chrom", pos_col: "pos"}, inplace=True)
if "chrom" not in df.columns:
    df["chrom"] = "NA"
if "pos" not in df.columns:
    df["pos"] = 0
df = df.dropna(subset=["pos", pval_col])
df = df[df[pval_col] > 0]
df["-log10p"] = -np.log10(df[pval_col])
df = df.sort_values(["chrom", "pos"])
x_pos = 0
for chrom in df["chrom"].unique():
    chrom_df = df[df["chrom"] == chrom]
    df.loc[chrom_df.index, "x"] = chrom_df["pos"] + x_pos
    x_pos += chrom_df["pos"].max() + 1e6  # 1Mb gap between chromosomes
plt.figure(figsize=(12, 6))
colors = ["skyblue", "navy"]
for i, chrom in enumerate(df["chrom"].unique()):
    chrom_df = df[df["chrom"] == chrom]
    plt.scatter(chrom_df["x"], chrom_df["-log10p"], color=colors[i % 2], s=10)
sig = df[df[pval_col] <= args.pvalue_threshold]
plt.scatter(sig["x"], sig["-log10p"], color="red", s=20, label="Significant")
plt.axhline(-np.log10(args.pvalue_threshold), color="red", linestyle="--", linewidth=1)
plt.xlabel("Chromosome")
plt.ylabel("-log10(p-value)")
plt.title("Manhattan Plot")
ticks = [df[df["chrom"] == chrom]["x"].median() for chrom in df["chrom"].unique()]
plt.xticks(ticks, df["chrom"].unique())
plt.tight_layout()
plt.savefig(args.plot_out, dpi=300)
plt.close()
top_hits = df.nsmallest(args.top_n, pval_col)
top_hits.to_csv(args.hits_out, sep="\t", index=False)

print(f"Manhattan plot saved to {args.plot_out}")
print(f"Top {args.top_n} hits saved to {args.hits_out}")
