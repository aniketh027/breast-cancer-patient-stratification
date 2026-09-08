# ============================================================================
# step13_compare.py
# Module Comparison: STRING vs Disease PPI Clusters
#
# Compares Leiden clustering from the STRING pipeline (step6) with the
# disease PPI pipeline (step11) using ARI, NMI, and Jaccard similarity.
#
# References:
#   [1] Hubert, L. & Arabie, P. (1985). Comparing partitions.
#       Journal of Classification, 2(1), 193–218.
#   [2] Vinh, N. X. et al. (2010). Information theoretic measures for
#       clusterings comparison. JMLR, 11, 2837–2854.
# ============================================================================

import pandas as pd
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
from itertools import product, combinations
import os

# ── Configuration ─────────────────────────────────────────────────────────────
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))

STRING_MODULES  = os.path.join(PROJECT_DIR, "results", "modules", "modules_BC_leiden.tsv")
DISEASE_MODULES = os.path.join(PROJECT_DIR, "results", "modules", "modules_BC_psicquic.tsv")
REPORT_DIR      = os.path.join(PROJECT_DIR, "results", "reports")

os.makedirs(REPORT_DIR, exist_ok=True)

print("=" * 60)
print("STEP 13: Module Comparison (STRING vs Disease PPI)")
print("=" * 60)

# ── 1. Load module assignments ───────────────────────────────────────────────
df_string = pd.read_csv(STRING_MODULES, sep="\t")
if "gene_symbol" in df_string.columns:
    df_string = df_string.rename(columns={"gene_symbol": "gene"})
if "cluster" in df_string.columns and "module_id" not in df_string.columns:
    df_string = df_string.rename(columns={"cluster": "module_id"})

df_disease = pd.read_csv(DISEASE_MODULES, sep="\t")
if "gene_symbol" in df_disease.columns:
    df_disease = df_disease.rename(columns={"gene_symbol": "gene"})
if "cluster" in df_disease.columns and "module_id" not in df_disease.columns:
    df_disease = df_disease.rename(columns={"cluster": "module_id"})

print(f"\nSTRING modules  : {df_string['module_id'].nunique()} modules, "
      f"{df_string['gene'].nunique()} unique genes")
print(f"Disease modules : {df_disease['module_id'].nunique()} modules, "
      f"{df_disease['gene'].nunique()} unique genes")

# ═══════════════════════════════════════════════════════════════════════════════
# 2. Jaccard module-pair similarity table
# ═══════════════════════════════════════════════════════════════════════════════
string_modules  = df_string.groupby("module_id")["gene"].apply(set).to_dict()
disease_modules = df_disease.groupby("module_id")["gene"].apply(set).to_dict()

jaccard_rows = []
for (s_mod, s_genes), (d_mod, d_genes) in product(
        string_modules.items(), disease_modules.items()):
    intersection = s_genes & d_genes
    union        = s_genes | d_genes
    jaccard = round(len(intersection) / len(union), 4) if union else 0.0
    jaccard_rows.append({
        "disease":         "BC",
        "string_module":   s_mod,
        "string_size":     len(s_genes),
        "disease_module":  d_mod,
        "disease_size":    len(d_genes),
        "overlap_size":    len(intersection),
        "jaccard":         jaccard,
        "overlap_genes":   ";".join(sorted(intersection)) if intersection else ""
    })

jaccard_df = pd.DataFrame(jaccard_rows).sort_values(
    ["string_module", "jaccard"], ascending=[True, False]
)
jaccard_df.to_csv(
    os.path.join(REPORT_DIR, "module_pair_similarity.tsv"),
    sep="\t", index=False
)

# Best match for each disease module
best_string_for_disease = (
    jaccard_df[jaccard_df["overlap_size"] > 0]
    .sort_values("jaccard", ascending=False)
    .drop_duplicates(subset="disease_module")
    .set_index("disease_module")["string_module"]
    .to_dict()
)

print(f"\nDisease PPI → STRING module mapping (by best Jaccard):")
for d_mod, s_mod in sorted(best_string_for_disease.items(),
                            key=lambda x: str(x[0])):
    print(f"  Disease {d_mod}  →  STRING {s_mod}")

# ═══════════════════════════════════════════════════════════════════════════════
# 3. Merge disease modules that map to the same STRING module
# ═══════════════════════════════════════════════════════════════════════════════
df_disease_merged = df_disease.copy()
df_disease_merged["module_id_merged"] = df_disease_merged["module_id"].map(
    lambda m: best_string_for_disease.get(m, m)
)

print(f"\nAfter merging disease modules to match STRING granularity:")
for merged_id, grp in df_disease_merged.groupby("module_id_merged"):
    orig_ids = sorted(grp["module_id"].unique())
    print(f"  Merged {orig_ids}  →  '{merged_id}'  ({len(grp)} genes)")

# ═══════════════════════════════════════════════════════════════════════════════
# 4. Compute ARI [1], NMI [2], Jaccard
# ═══════════════════════════════════════════════════════════════════════════════
merged_raw = pd.merge(
    df_string[["gene", "module_id"]],
    df_disease[["gene", "module_id"]],
    on="gene", suffixes=("_string", "_disease")
)

merged_fixed = pd.merge(
    df_string[["gene", "module_id"]],
    df_disease_merged[["gene", "module_id_merged"]],
    on="gene"
)

print(f"\nShared universe (raw)   : {len(merged_raw)} genes")
print(f"Shared universe (merged): {len(merged_fixed)} genes")

report_rows = []

# Jaccard clustering similarity
def jaccard_clustering(labels1, labels2):
    n = len(labels1)
    both = 0
    either = 0
    for i, j in combinations(range(n), 2):
        same_1 = labels1[i] == labels1[j]
        same_2 = labels2[i] == labels2[j]
        if same_1 and same_2:
            both += 1
        if same_1 or same_2:
            either += 1
    return both / either if either > 0 else 0.0

if len(merged_raw) >= 2:
    nmi_raw = normalized_mutual_info_score(
        merged_raw["module_id_string"], merged_raw["module_id_disease"])
    ari_raw = adjusted_rand_score(
        merged_raw["module_id_string"], merged_raw["module_id_disease"])
    jac_raw = jaccard_clustering(
        merged_raw["module_id_string"].tolist(),
        merged_raw["module_id_disease"].tolist()
    )
    report_rows += [
        {"disease": "BC", "viewA": "string_leiden", "viewB": "disease_leiden",
         "comparison": "raw", "metric": "NMI",
         "universe_size": len(merged_raw), "value": round(nmi_raw, 4)},
        {"disease": "BC", "viewA": "string_leiden", "viewB": "disease_leiden",
         "comparison": "raw", "metric": "ARI",
         "universe_size": len(merged_raw), "value": round(ari_raw, 4)},
        {"disease": "BC", "viewA": "string_leiden", "viewB": "disease_leiden",
         "comparison": "raw", "metric": "Jaccard",
         "universe_size": len(merged_raw), "value": round(jac_raw, 4)},
    ]

if len(merged_fixed) >= 2:
    nmi_merged = normalized_mutual_info_score(
        merged_fixed["module_id"], merged_fixed["module_id_merged"])
    ari_merged = adjusted_rand_score(
        merged_fixed["module_id"], merged_fixed["module_id_merged"])
    report_rows += [
        {"disease": "BC", "viewA": "string_leiden", "viewB": "disease_merged",
         "comparison": "granularity_matched", "metric": "NMI",
         "universe_size": len(merged_fixed), "value": round(nmi_merged, 4)},
        {"disease": "BC", "viewA": "string_leiden", "viewB": "disease_merged",
         "comparison": "granularity_matched", "metric": "ARI",
         "universe_size": len(merged_fixed), "value": round(ari_merged, 4)},
    ]

report_df = pd.DataFrame(report_rows)
report_df.to_csv(
    os.path.join(REPORT_DIR, "module_comparison.tsv"),
    sep="\t", index=False
)

# ── 5. Print summary ─────────────────────────────────────────────────────────
print(f"\n{'─' * 50}")
print("Comparison Metrics:")
print(f"{'─' * 50}")
if not report_df.empty:
    print(report_df.to_string(index=False))
else:
    print("  No overlapping genes found between STRING and Disease PPI modules.")

if report_rows:
    ari_r = next((r["value"] for r in report_rows
                  if r["metric"] == "ARI" and r["comparison"] == "raw"), None)
    ari_m = next((r["value"] for r in report_rows
                  if r["metric"] == "ARI"
                  and r["comparison"] == "granularity_matched"), None)
    if ari_r is not None and ari_m is not None:
        print(f"\nARI improvement from granularity matching: "
              f"{ari_r:.4f} → {ari_m:.4f} (+{ari_m - ari_r:.4f})")

# ── 6. Best match summary ────────────────────────────────────────────────────
print(f"\n{'─' * 50}")
print("Best disease PPI match for each STRING module:")
print(f"{'─' * 50}")
best_matches = (
    jaccard_df[jaccard_df["overlap_size"] > 0]
    .sort_values("jaccard", ascending=False)
    .drop_duplicates(subset="string_module")
    .sort_values("string_module")
)
for _, row in best_matches.iterrows():
    print(f"  STRING {row['string_module']} (n={row['string_size']}) "
          f"<-> Disease {row['disease_module']} (n={row['disease_size']}) "
          f"| overlap={row['overlap_size']}  Jaccard={row['jaccard']}")
    if row["overlap_genes"]:
        print(f"    shared: {row['overlap_genes']}")

print(f"\n{'─' * 50}")
print(f"Saved: {os.path.join(REPORT_DIR, 'module_comparison.tsv')}")
print(f"Saved: {os.path.join(REPORT_DIR, 'module_pair_similarity.tsv')}")
print(f"\n✅ Step 13 complete — Module comparison done")
