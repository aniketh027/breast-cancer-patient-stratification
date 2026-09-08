# ============================================================================
# step8_disease_ppi.py
# Disease-Specific PPI Network via PSICQUIC + BioGRID
#
# Builds a disease-specific PPI by querying IntAct via PSICQUIC using
# the disease term "breast cancer" (NOT by filtering around seed genes).
# Also integrates BioGRID physical interactions for genes found in the
# IntAct disease query.
#
# This produces a disease-contextualized network independent of seed genes,
# allowing a fair comparison with the STRING-based pipeline.
#
# References:
#   [1] del-Toro, N. et al. (2013). A new reference implementation of the
#       PSICQUIC web service. Nucleic Acids Research, 41(W1), W601–W606.
#   [2] Oughtred, R. et al. (2021). The BioGRID database. Nucleic Acids Res.
# ============================================================================

import pandas as pd
import numpy as np
import requests
import zipfile
import subprocess
import re
import os
import time
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────────
PROJECT_DIR    = os.path.dirname(os.path.abspath(__file__))
DISEASE_ABBR   = "BC"
DISEASE_TERMS  = ["breast cancer", "breast neoplasm", "breast carcinoma"]

SEEDS_FILE     = os.path.join(PROJECT_DIR, "results", "genes", "seeds_BC.tsv")
BIOGRID_ZIP    = os.path.join(PROJECT_DIR, "resources", "biogrid",
                               "BIOGRID-MV-Physical-LATEST.tab3.zip")
OUT_FILE       = os.path.join(PROJECT_DIR, "results", "networks",
                               f"ppi_disease_{DISEASE_ABBR}.tsv")
BIOGRID_URL    = ("https://downloads.thebiogrid.org/Download/BioGRID/"
                  "Latest-Release/BIOGRID-MV-Physical-LATEST.tab3.zip")

MI_SCORE_CUTOFF = 0.4   # IntAct MI score threshold
PSICQUIC_BATCH  = 2000  # records per PSICQUIC page
MAX_RECORDS     = 100000  # cap to keep download time reasonable

os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)

print("=" * 60)
print(f"STEP 8: Disease-Specific PPI Network (PSICQUIC + BioGRID)")
print("=" * 60)

# ── Load seeds (for reporting only — NOT used for filtering) ──────────────────
seeds_df   = pd.read_csv(SEEDS_FILE, sep="\t")
seed_set   = set(seeds_df["gene_symbol"].str.strip().str.upper())
print(f"\nSeed genes (for reference): {len(seed_set)}")

# ══════════════════════════════════════════════════════════════════════════════
# SOURCE 1: IntAct via PSICQUIC [1]
# Query by disease term, NOT by seed genes
# ══════════════════════════════════════════════════════════════════════════════
PSICQUIC_BASE = ("https://www.ebi.ac.uk/Tools/webservices/psicquic/"
                 "intact/webservices/current/search/query/")

GENE_PAT = re.compile(r'(?:uniprotkb|hgnc):([A-Z][A-Z0-9\-]+)\(gene name\)', re.I)

def query_psicquic(disease_term, max_records=MAX_RECORDS):
    """Query IntAct PSICQUIC for human-human interactions with disease term."""
    query = f'species:9606 AND species:9606 AND disease:"{disease_term}"'
    encoded = requests.utils.quote(query)

    print(f"\n  Querying PSICQUIC: {disease_term}")

    # Get count first
    count_url = f"{PSICQUIC_BASE}{encoded}?format=count"
    try:
        resp = requests.get(count_url, timeout=30)
        total = int(resp.text.strip())
        print(f"  Total results: {total}")
    except Exception as e:
        print(f"  ❌ Count failed: {e}")
        return []

    # Download in batches
    records = []
    n_download = min(total, max_records)
    print(f"  Downloading up to {n_download} records...")

    for start in range(0, n_download, PSICQUIC_BATCH):
        url = (f"{PSICQUIC_BASE}{encoded}"
               f"?format=tab25&firstResult={start}&maxResults={PSICQUIC_BATCH}")
        try:
            resp = requests.get(url, timeout=60)
            if resp.status_code != 200:
                print(f"    ⚠ HTTP {resp.status_code} at offset {start}")
                continue

            lines = resp.text.strip().split("\n")
            batch_kept = 0

            for line in lines:
                if not line.strip():
                    continue
                fields = line.split("\t")
                if len(fields) < 15:
                    continue

                # Parse taxids (columns 9,10 in MITAB 2.5)
                tax_a, tax_b = fields[9], fields[10]
                if "9606" not in tax_a or "9606" not in tax_b:
                    continue   # skip non-human

                # Parse gene names from aliases (columns 4,5)
                alias_a, alias_b = fields[4], fields[5]
                genes_a = GENE_PAT.findall(alias_a)
                genes_b = GENE_PAT.findall(alias_b)

                if not genes_a or not genes_b:
                    continue

                gene_a = genes_a[0].upper()
                gene_b = genes_b[0].upper()

                if gene_a == gene_b:
                    continue   # skip self-loops

                # Parse MI score (column 14)
                conf = fields[14] if len(fields) > 14 else ""
                mi_match = re.search(r'intact-miscore:([\d.]+)', conf)
                weight = float(mi_match.group(1)) if mi_match else 0.0

                if weight < MI_SCORE_CUTOFF:
                    continue

                records.append({
                    "geneA": gene_a,
                    "geneB": gene_b,
                    "weight": weight
                })
                batch_kept += 1

            if (start // PSICQUIC_BATCH) % 10 == 0:
                print(f"    Batch {start}-{start+PSICQUIC_BATCH}: "
                      f"{len(lines)} raw → {batch_kept} kept")

            time.sleep(0.2)   # be nice to the API

        except Exception as e:
            print(f"    ⚠ Error at offset {start}: {e}")
            continue

    print(f"  ✅ IntAct PSICQUIC: {len(records)} interactions kept")
    return records


# Query with multiple disease terms and combine
all_intact = []
for term in DISEASE_TERMS:
    results = query_psicquic(term)
    all_intact.extend(results)

if all_intact:
    df_intact = pd.DataFrame(all_intact)
    df_intact["source"] = "IntAct"
    # Deduplicate within IntAct (same edge from different terms)
    df_intact["pair"] = df_intact.apply(
        lambda r: tuple(sorted([r["geneA"], r["geneB"]])), axis=1)
    df_intact = (df_intact.groupby("pair")
                 .agg(geneA=("geneA", "first"),
                      geneB=("geneB", "first"),
                      weight=("weight", "max"),
                      source=("source", "first"))
                 .reset_index(drop=True))
    print(f"\n  IntAct deduplicated: {len(df_intact)} unique edges")
else:
    df_intact = pd.DataFrame(columns=["geneA", "geneB", "weight", "source"])
    print("\n  ⚠ No IntAct results")

# ══════════════════════════════════════════════════════════════════════════════
# SOURCE 2: BioGRID [2]
# Use ALL human physical interactions (BioGRID doesn't support disease queries)
# ══════════════════════════════════════════════════════════════════════════════
def download_biogrid():
    """Download BioGRID if not present."""
    if os.path.exists(BIOGRID_ZIP):
        try:
            zipfile.ZipFile(BIOGRID_ZIP).testzip()
            print(f"  ✅ BioGRID exists: {os.path.basename(BIOGRID_ZIP)}")
            return
        except:
            os.remove(BIOGRID_ZIP)

    print(f"  ⬇  Downloading BioGRID...")
    os.makedirs(os.path.dirname(BIOGRID_ZIP), exist_ok=True)
    subprocess.run(["curl", "-L", "-o", BIOGRID_ZIP, BIOGRID_URL],
                   check=True, capture_output=True)
    print(f"  ✅ Downloaded ({os.path.getsize(BIOGRID_ZIP)/1e6:.1f} MB)")


def parse_biogrid(zip_path, disease_genes):
    """Parse BioGRID for human interactions involving disease-relevant genes."""
    print(f"\nParsing BioGRID (filtered by {len(disease_genes)} disease genes)...")
    results = []

    with zipfile.ZipFile(zip_path) as z:
        tab3_file = next(
            (n for n in z.namelist() if n.endswith(".tab3.txt")), None)
        if not tab3_file:
            print("  ❌ No TAB3 file found")
            return pd.DataFrame(columns=["geneA", "geneB", "weight"])

        with z.open(tab3_file) as f:
            chunks = pd.read_csv(f, sep="\t", chunksize=50000,
                                 low_memory=False, header=0)

            for chunk_num, chunk in enumerate(chunks):
                chunk.columns = [c.strip().lstrip("#").strip()
                                 for c in chunk.columns]

                sym_a = next((c for c in chunk.columns
                              if "official symbol" in c.lower()
                              and "interactor a" in c.lower()), None)
                sym_b = next((c for c in chunk.columns
                              if "official symbol" in c.lower()
                              and "interactor b" in c.lower()), None)
                org_a = next((c for c in chunk.columns
                              if "organism" in c.lower()
                              and "interactor a" in c.lower()), None)
                org_b = next((c for c in chunk.columns
                              if "organism" in c.lower()
                              and "interactor b" in c.lower()), None)

                if not sym_a or not sym_b:
                    continue

                # Human filter
                if org_a and org_b:
                    human = (
                        chunk[org_a].astype(str).str.contains(
                            "9606|Homo sapiens", na=False, case=False) &
                        chunk[org_b].astype(str).str.contains(
                            "9606|Homo sapiens", na=False, case=False)
                    )
                    chunk = chunk[human].copy()

                if chunk.empty:
                    continue

                chunk["geneA"] = chunk[sym_a].astype(str).str.strip().str.upper()
                chunk["geneB"] = chunk[sym_b].astype(str).str.strip().str.upper()
                chunk["weight"] = 1.0   # BioGRID: uniform weight

                chunk = chunk[~chunk["geneA"].isin(["-", "nan", ""])]
                chunk = chunk[~chunk["geneB"].isin(["-", "nan", ""])]
                chunk = chunk[chunk["geneA"] != chunk["geneB"]]

                # Disease gene filter — BOTH genes must be in disease set
                disease_mask = (chunk["geneA"].isin(disease_genes) &
                                chunk["geneB"].isin(disease_genes))
                chunk = chunk[disease_mask].copy()

                if chunk.empty:
                    continue

                results.append(chunk[["geneA", "geneB", "weight"]])

    if not results:
        return pd.DataFrame(columns=["geneA", "geneB", "weight"])

    df = pd.concat(results, ignore_index=True)
    df["source"] = "BioGRID"
    print(f"  ✅ BioGRID disease-filtered edges: {len(df)}")
    return df


# Extract disease gene set from IntAct PSICQUIC results
intact_genes = set()
if len(df_intact) > 0:
    intact_genes = set(df_intact["geneA"]) | set(df_intact["geneB"])
print(f"\n  Disease gene set from IntAct: {len(intact_genes)} genes")

print(f"\n{'='*50}")
print("Checking BioGRID data...")
download_biogrid()
df_biogrid = parse_biogrid(BIOGRID_ZIP, intact_genes)

# ══════════════════════════════════════════════════════════════════════════════
# MERGE IntAct + BioGRID
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*50}")
print(f"IntAct edges  : {len(df_intact)}")
print(f"BioGRID edges : {len(df_biogrid)}")

df_all = pd.concat([df_intact, df_biogrid], ignore_index=True)

# Normalize to undirected pairs
df_all["pair"] = df_all.apply(
    lambda r: tuple(sorted([r["geneA"], r["geneB"]])), axis=1)

# For duplicates — keep max weight
df_merged = (df_all.groupby("pair")
             .agg(weight=("weight", "max"),
                  sources=("source", lambda x: "+".join(sorted(set(x)))))
             .reset_index()
             .assign(
                 geneA=lambda x: x["pair"].apply(lambda p: p[0]),
                 geneB=lambda x: x["pair"].apply(lambda p: p[1])
             )[["geneA", "geneB", "weight", "sources"]])

# Rename columns for pipeline compatibility
df_merged = df_merged.rename(columns={"geneA": "nodeA", "geneB": "nodeB"})

# Summary
all_genes     = set(df_merged["nodeA"]) | set(df_merged["nodeB"])
seeds_covered = seed_set & all_genes
both_sources  = df_merged[df_merged["sources"].str.contains(r"\+")]

print(f"\n✅ Total merged edges     : {len(df_merged)}")
print(f"✅ Unique genes           : {len(all_genes)}")
print(f"✅ Seeds in network       : {len(seeds_covered)}/{len(seed_set)}")
print(f"✅ Edges in BOTH sources  : {len(both_sources)}  ← high confidence")

# Save
Path(OUT_FILE).parent.mkdir(parents=True, exist_ok=True)
df_merged.to_csv(OUT_FILE, sep="\t", index=False)
print(f"\n✅ Saved → {OUT_FILE}")

