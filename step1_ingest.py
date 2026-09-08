import pandas as pd
import os

os.makedirs("results/genes", exist_ok=True)

raw_df = pd.read_csv("resources/raw/gwas/BC_reported_genes.tsv", sep="\t")

raw_df["source"] = "BC_latest_reported_genes"
raw_df["evidence"] = "reported_gwas_gene"

raw_df.to_csv("results/genes/seeds_BC.tsv", sep="\t", index=False)
