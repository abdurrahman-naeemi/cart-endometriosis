"""
Tabula Sapiens — Double Stratified Downsampling (RAM Optimised)
----------------------------------------------------------------
Key RAM optimisation: never holds all chunks in memory at once.
Writes expression data directly to output h5ad incrementally
using h5py, so peak RAM usage stays low throughout.

Usage:
    python downsample_tabula_sapiens.py

Output:
    tabula_sapiens_double_stratified.h5ad
"""

import h5py
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
import os
import gc

# ── CONFIG ────────────────────────────────────────────────────────────────────
INPUT_PATH  = "./tabula_sapiens.h5ad"
OUTPUT_PATH = "./tabula_sapiens_double_stratified.h5ad"
MAX_CELLS_PER_COMBO = 5000
DONOR_COL   = "donor_id"
TISSUE_COL  = "tissue"
SEED        = 42
CHUNK_SIZE  = 5000       # rows to read at once from input
WRITE_CHUNK = 5000      # rows to write at once to output

np.random.seed(SEED)

# ── STEP 1: LOAD OBS ONLY ─────────────────────────────────────────────────────
print("Reading cell metadata (obs only)...")
adata_backed = sc.read_h5ad(INPUT_PATH, backed="r")
obs = adata_backed.obs.copy()
n_total = len(obs)
adata_backed.file.close()
del adata_backed
gc.collect()

print(f"  Total cells:    {n_total:,}")
print(f"  Unique donors:  {obs[DONOR_COL].nunique()}")
print(f"  Unique tissues: {obs[TISSUE_COL].nunique()}")

# ── STEP 2: DOUBLE STRATIFIED SAMPLING ────────────────────────────────────────
print(f"\nSampling up to {MAX_CELLS_PER_COMBO} cells per donor x tissue combo...")

sampled_indices = []
combo_stats = []

for (donor, tissue), group in obs.groupby([DONOR_COL, TISSUE_COL], observed=True):
    n_available = len(group)
    n_sample    = min(n_available, MAX_CELLS_PER_COMBO)
    sampled     = group.sample(n=n_sample, replace=False, random_state=SEED)
    sampled_indices.extend(sampled.index.tolist())
    combo_stats.append({
        "donor": donor, "tissue": tissue,
        "available": n_available, "sampled": n_sample
    })

combo_df = pd.DataFrame(combo_stats)
total_sampled = len(sampled_indices)
print(f"  Total cells to extract: {total_sampled:,}")

combo_df.to_csv("sampling_stats.csv", index=False)

# Get sorted integer positions
all_barcodes   = obs.index.tolist()
barcode_to_pos = {bc: i for i, bc in enumerate(all_barcodes)}
sampled_pos    = np.array(sorted([barcode_to_pos[bc] for bc in sampled_indices]))

# Save obs for sampled cells now (small, fine in RAM)
obs_final = obs.iloc[sampled_pos].copy()

del obs, all_barcodes, barcode_to_pos, sampled_indices
gc.collect()

# ── STEP 3: STREAM + WRITE INCREMENTALLY ─────────────────────────────────────
# Instead of collecting all chunks then vstacking (kills RAM),
# we write each batch directly to the output h5py file as we go.
# Final file is assembled on disk, not in RAM.

print(f"\nStreaming input and writing output incrementally...")
print(f"  (peak RAM stays low — no vstack)")

# First pass: collect all sampled rows as list of small sparse matrices
# but flush to disk every WRITE_CHUNK rows
tmp_path = OUTPUT_PATH + ".tmp.h5"

with h5py.File(INPUT_PATH, "r") as f_in:
    x_group = f_in["X"]
    is_sparse = isinstance(x_group, h5py.Group)

    if not is_sparse:
        raise ValueError("Expected sparse CSR matrix in input h5ad")

    data_ds    = x_group["data"]
    indices_ds = x_group["indices"]
    indptr_ds  = x_group["indptr"]
    n_rows     = len(indptr_ds) - 1
    n_cols     = len(f_in["var"]["_index"])
    print(f"  Input matrix: {n_rows:,} x {n_cols:,}")

    # Pre-allocate output arrays — write to temp h5 file incrementally
    with h5py.File(tmp_path, "w") as f_out:
        # We'll write CSR: data, indices, indptr
        # Use resizable datasets
        ds_data    = f_out.create_dataset("data",    shape=(0,), maxshape=(None,),
                                          dtype=np.float32, chunks=(1000000,))
        ds_indices = f_out.create_dataset("indices", shape=(0,), maxshape=(None,),
                                          dtype=np.int32,   chunks=(1000000,))
        ds_indptr  = f_out.create_dataset("indptr",  shape=(1,), maxshape=(None,),
                                          dtype=np.int64,   chunks=(1000,))
        ds_indptr[0] = 0

        current_nnz   = 0
        rows_written  = 0
        chunk_start   = 0
        chunk_count   = 0

        while chunk_start < n_rows:
            chunk_end = min(chunk_start + CHUNK_SIZE, n_rows)

            in_chunk = sampled_pos[
                (sampled_pos >= chunk_start) &
                (sampled_pos < chunk_end)
            ]

            if len(in_chunk) > 0:
                row_start = int(indptr_ds[chunk_start])
                row_end   = int(indptr_ds[chunk_end])

                chunk_data    = data_ds[row_start:row_end]
                chunk_indices = indices_ds[row_start:row_end]
                chunk_indptr  = indptr_ds[chunk_start:chunk_end+1] - indptr_ds[chunk_start]

                chunk_matrix = sp.csr_matrix(
                    (chunk_data, chunk_indices, chunk_indptr),
                    shape=(chunk_end - chunk_start, n_cols)
                )

                local_pos   = in_chunk - chunk_start
                sampled_mat = chunk_matrix[local_pos, :]

                # Write to output file immediately
                new_data    = sampled_mat.data.astype(np.float32)
                new_indices = sampled_mat.indices.astype(np.int32)
                new_indptr  = sampled_mat.indptr[1:].astype(np.int64) + current_nnz

                new_nnz = len(new_data)
                ds_data.resize((current_nnz + new_nnz,))
                ds_indices.resize((current_nnz + new_nnz,))
                ds_data[current_nnz:]    = new_data
                ds_indices[current_nnz:] = new_indices

                n_new_rows = sampled_mat.shape[0]
                ds_indptr.resize((rows_written + n_new_rows + 1,))
                ds_indptr[rows_written + 1: rows_written + n_new_rows + 1] = new_indptr

                current_nnz  += new_nnz
                rows_written += n_new_rows

                # Free immediately
                del chunk_matrix, sampled_mat, chunk_data, chunk_indices, chunk_indptr
                del new_data, new_indices, new_indptr
                gc.collect()

                chunk_count += 1
                if chunk_count % 20 == 0:
                    print(f"  Chunk {chunk_count} | rows written: {rows_written:,} / {total_sampled:,}")

            chunk_start = chunk_end

        f_out.attrs["shape"] = [rows_written, n_cols]
        print(f"  Done streaming. Total rows: {rows_written:,}, nnz: {current_nnz:,}")

# ── STEP 4: ASSEMBLE FINAL H5AD FROM TEMP FILE ───────────────────────────────
print("\nAssembling final h5ad from temp file...")
print("  (reading back in chunks — stays under RAM limit)")

# Read back in chunks and build final AnnData
chunk_matrices = []
with h5py.File(tmp_path, "r") as f_tmp:
    data_arr    = f_tmp["data"][:]
    indices_arr = f_tmp["indices"][:]
    indptr_arr  = f_tmp["indptr"][:]
    shape       = tuple(f_tmp.attrs["shape"])

X_final = sp.csr_matrix((data_arr, indices_arr, indptr_arr), shape=shape)
del data_arr, indices_arr, indptr_arr
gc.collect()
print(f"  Final matrix: {X_final.shape}")

# Get var metadata
adata_meta = sc.read_h5ad(INPUT_PATH, backed="r")
var_final  = adata_meta.var.copy()
adata_meta.file.close()
del adata_meta
gc.collect()

adata_out = ad.AnnData(X=X_final, obs=obs_final, var=var_final)
del X_final
gc.collect()

print(f"\nFinal AnnData: {adata_out.n_obs:,} cells x {adata_out.n_vars:,} genes")
print(f"\nDonor distribution:")
print(adata_out.obs[DONOR_COL].value_counts().to_string())

# ── STEP 5: SAVE ──────────────────────────────────────────────────────────────
print(f"\nSaving to {OUTPUT_PATH}...")
adata_out.write_h5ad(OUTPUT_PATH, compression="gzip")
size_gb = os.path.getsize(OUTPUT_PATH) / 1e9
print(f"Done. File size: {size_gb:.2f} GB")

# Clean up temp file
os.remove(tmp_path)
print(f"Cleaned up temp file.")
print(f"\nUpload {OUTPUT_PATH} to OneDrive and share with team.")