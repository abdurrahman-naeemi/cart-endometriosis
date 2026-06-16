"""
Surfaceome x scRNA-seq intersection
------------------------------------
Takes the Bausch-Fluck 2886 surfaceome CSV and intersects it with
the genes present in your combined h5ad to give a shortlist of
cell surface marker candidates.

Usage:
    python intersect_markers.py

Outputs:
    surfaceome_in_atlas.csv   — full rows from surfaceome CSV for matched genes
    marker_shortlist.txt      — clean gene list, one per line
"""

import pandas as pd
import scanpy as sc
import os

# ── CONFIG ────────────────────────────────────────────────────────────────────
H5AD_PATH       = "./all_samples_combined_w_meta_QC_scvi_umap.h5ad"
SURFACEOME_PATH = "./surfaceome.csv"
OUT_DIR         = "./marker_outputs"

os.makedirs(OUT_DIR, exist_ok=True)

# ── LOAD SURFACEOME ───────────────────────────────────────────────────────────
print("Loading surfaceome CSV...")
surf = pd.read_csv(SURFACEOME_PATH, skiprows=1)   # skip the title row
surf.columns = surf.columns.str.strip()

# Gene names are in the 'UniProt gene' column
gene_col = "UniProt gene"
surf = surf.dropna(subset=[gene_col])
surf[gene_col] = surf[gene_col].str.strip().str.upper()

surfaceome_genes = set(surf[gene_col].tolist())
print(f"  Surfaceome: {len(surfaceome_genes)} genes")

# ── LOAD H5AD ─────────────────────────────────────────────────────────────────
print("Loading h5ad...")
adata = sc.read_h5ad(H5AD_PATH)
atlas_genes = set(adata.var_names.str.upper().tolist())
print(f"  Atlas: {adata.n_obs} cells x {adata.n_vars} genes")

# ── INTERSECT ─────────────────────────────────────────────────────────────────
overlap = surfaceome_genes & atlas_genes
print(f"\n✓ Overlap: {len(overlap)} surface markers found in atlas")
print(f"  ({len(surfaceome_genes) - len(overlap)} surfaceome genes not in atlas — likely not expressed / filtered out)")

# ── BUILD OUTPUT TABLE ────────────────────────────────────────────────────────
# Pull full surfaceome metadata for matched genes
matched = surf[surf[gene_col].isin(overlap)].copy()

# Useful columns to keep
keep_cols = [
    "UniProt gene",
    "UniProt description",
    "Surfaceome Label",
    "Surfaceome Label Source",
    "CD number",
    "CSPA category",          # mass spec validation confidence
    "TM domains",
    "topology",
    "HPA antibody",           # whether an antibody already exists (useful for CAR design)
    "DrugBank approved drug IDs",
    "UniProt subcellular",
]

keep_cols = [c for c in keep_cols if c in matched.columns]
matched_clean = matched[keep_cols].sort_values("UniProt gene")

# ── CSPA FILTER (optional but recommended) ───────────────────────────────────
# CSPA = Cell Surface Protein Atlas — mass spec validated
# Category 1 = high confidence, Category 2 = medium
# Filter to CSPA validated only for highest quality shortlist
cspa_validated = matched_clean[
    matched_clean["CSPA category"].astype(str).str.contains("1|2", na=False)
] if "CSPA category" in matched_clean.columns else matched_clean

print(f"\n  CSPA mass-spec validated subset: {len(cspa_validated)} genes")
print(f"  (these have direct experimental evidence of surface expression)")

# ── SAVE ──────────────────────────────────────────────────────────────────────
# Full intersection
full_path = os.path.join(OUT_DIR, "surfaceome_in_atlas.csv")
matched_clean.to_csv(full_path, index=False)
print(f"\nSaved full intersection → {full_path}")

# CSPA validated only
cspa_path = os.path.join(OUT_DIR, "surfaceome_in_atlas_CSPA_validated.csv")
cspa_validated.to_csv(cspa_path, index=False)
print(f"Saved CSPA-validated subset → {cspa_path}")

# Plain gene list
genes_path = os.path.join(OUT_DIR, "marker_shortlist.txt")
with open(genes_path, "w") as f:
    f.write("\n".join(sorted(overlap)))
print(f"Saved plain gene list → {genes_path}")

# ── PREVIEW ───────────────────────────────────────────────────────────────────
print("\n── Top 20 CSPA-validated candidates ──────────────────────────────")
preview_cols = ["UniProt gene", "UniProt description", "CSPA category", "CD number", "HPA antibody"]
preview_cols = [c for c in preview_cols if c in cspa_validated.columns]
print(cspa_validated[preview_cols].head(20).to_string(index=False))

print("\nDone. Next step: load surfaceome_in_atlas_CSPA_validated.csv into")
print("your scVI pipeline and score each gene for ectopic vs eutopic expression.")
