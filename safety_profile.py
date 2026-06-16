"""
Human Safety Profile v2 — GAS1 + BST2
----------------------------------------
Upgrades from v1:
  1. Threshold sweep (0.2, 0.5, 1.0, 2.0, 5.0) instead of single threshold
  2. Transfer curves per organ — co-expression rate vs threshold
  3. Patient-stratified rates — compute per donor, then average across donors

Usage:
    python safety_profile_v2.py

Output:
    safety_outputs_v2/transfer_curves.png      — co-expression vs threshold per tissue
    safety_outputs_v2/safety_report_t1.csv     — full report at threshold >1
    safety_outputs_v2/patient_stratified.csv   — patient-stratified rates at all thresholds
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.sparse as sp
import os

# ── CONFIG ────────────────────────────────────────────────────────────────────
TABULA_PATH   = "./tabula_sapiens_10pct.h5ad"
OUT_DIR       = "./safety_outputs_v2"

# Ensembl IDs
TARGET_GENES  = ["ENSG00000180447", "ENSG00000130303"]
TARGET_LABELS = ["GAS1",            "BST2"]

# Threshold sweep — Haziq wants >1 as primary, but sweep all of these
THRESHOLDS = [0.2, 0.5, 1.0, 2.0, 5.0]
PRIMARY_THRESHOLD = 1.0   # the one Haziq asked for

# Critical organs
CRITICAL_ORGANS = [
    "heart", "lung", "brain", "kidney", "liver",
    "spinal cord", "cardiac muscle", "cerebellum", "cerebral cortex"
]

os.makedirs(OUT_DIR, exist_ok=True)

# ── LOAD ──────────────────────────────────────────────────────────────────────
print("Loading Tabula Sapiens...")
ts = sc.read_h5ad(TABULA_PATH)
print(f"  {ts.n_obs:,} cells x {ts.n_vars:,} genes")
print(f"  Tissues: {ts.obs['tissue'].nunique()}")
print(f"  Donors:  {ts.obs['donor_id'].nunique()}")

# ── FIND GENES ────────────────────────────────────────────────────────────────
print(f"\nLooking for target genes...")
found = []
for ensembl, label in zip(TARGET_GENES, TARGET_LABELS):
    if ensembl in ts.var_names:
        found.append((ensembl, label))
        print(f"  ✓ {label} ({ensembl})")
    else:
        print(f"  ✗ {label} ({ensembl}) NOT FOUND")

if not found:
    raise ValueError("No target genes found.")

labels = [l for _, l in found]

# ── NORMALISE ─────────────────────────────────────────────────────────────────
print("\nChecking normalisation...")
sample = ts.X[:100, :]
max_val = sample.max() if not sp.issparse(sample) else sample.max()
if max_val > 100:
    print(f"  Normalising (max={max_val:.0f})...")
    sc.pp.normalize_total(ts, target_sum=1e4)
    sc.pp.log1p(ts)
else:
    print(f"  Already normalised (max={max_val:.2f})")

# ── EXTRACT EXPRESSION ────────────────────────────────────────────────────────
print(f"\nExtracting expression for {labels}...")
expr = {}
for ensembl, label in found:
    idx = ts.var_names.get_loc(ensembl)
    col = ts.X[:, idx]
    if sp.issparse(col):
        col = np.asarray(col.todense()).flatten()
    expr[label] = col

tissues   = ts.obs["tissue"].values
donors    = ts.obs["donor_id"].values
tissue_list = sorted(ts.obs["tissue"].unique())
donor_list  = sorted(ts.obs["donor_id"].unique())

# ── PATIENT-STRATIFIED RATES ──────────────────────────────────────────────────
# For each tissue + threshold:
#   1. Compute co-expression rate within each donor that has cells in that tissue
#   2. Average across donors (each donor weighted equally)
# This prevents one outlier donor from dominating

print(f"\nComputing patient-stratified rates across {len(THRESHOLDS)} thresholds...")
print(f"  Tissues: {len(tissue_list)}  x  Thresholds: {len(THRESHOLDS)}")

# Store results: dict[tissue][threshold] = {mean, std, n_donors, per_donor_rates}
strat_results = {}

for tissue in tissue_list:
    strat_results[tissue] = {}
    tissue_mask = tissues == tissue

    for threshold in THRESHOLDS:
        donor_coexpr_rates = []

        for donor in donor_list:
            donor_mask = donors == donor
            cell_mask  = tissue_mask & donor_mask
            n_cells    = cell_mask.sum()

            if n_cells < 5:   # skip donors with too few cells in this tissue
                continue

            # Co-expression in this donor x tissue
            activated = {}
            for label in labels:
                activated[label] = expr[label][cell_mask] > threshold

            coexpr = activated[labels[0]]
            for label in labels[1:]:
                coexpr = coexpr & activated[label]

            donor_coexpr_rates.append(float(coexpr.mean()) * 100)

        if donor_coexpr_rates:
            strat_results[tissue][threshold] = {
                "mean":     np.mean(donor_coexpr_rates),
                "std":      np.std(donor_coexpr_rates),
                "n_donors": len(donor_coexpr_rates),
                "rates":    donor_coexpr_rates,
            }
        else:
            strat_results[tissue][threshold] = {
                "mean": 0.0, "std": 0.0, "n_donors": 0, "rates": []
            }

# ── BUILD REPORT AT PRIMARY THRESHOLD ────────────────────────────────────────
print(f"\nBuilding report at primary threshold >{PRIMARY_THRESHOLD}...")
report_rows = []
for tissue in tissue_list:
    r = strat_results[tissue].get(PRIMARY_THRESHOLD, {})
    is_critical = any(o.lower() in tissue.lower() for o in CRITICAL_ORGANS)
    mean_rate   = r.get("mean", 0.0)
    report_rows.append({
        "tissue":            tissue,
        "n_donors":          r.get("n_donors", 0),
        "coexpr_mean_pct":   round(mean_rate, 3),
        "coexpr_std_pct":    round(r.get("std", 0.0), 3),
        "is_critical_organ": is_critical,
        "safety_flag":       "CRITICAL" if (is_critical and mean_rate > 5) else
                             "WARNING"  if mean_rate > 10 else
                             "SAFE"
    })

report_df = pd.DataFrame(report_rows).sort_values("coexpr_mean_pct", ascending=False)
report_df.to_csv(os.path.join(OUT_DIR, "safety_report_t1.csv"), index=False)

# ── BUILD TRANSFER CURVE TABLE ────────────────────────────────────────────────
print("Building transfer curve table...")
curve_rows = []
for tissue in tissue_list:
    row = {"tissue": tissue}
    for t in THRESHOLDS:
        r = strat_results[tissue].get(t, {})
        row[f"t{t}"] = round(r.get("mean", 0.0), 3)
    curve_rows.append(row)

curve_df = pd.DataFrame(curve_rows)
curve_df.to_csv(os.path.join(OUT_DIR, "patient_stratified.csv"), index=False)

# ── PRINT REPORT ──────────────────────────────────────────────────────────────
print("\n" + "="*70)
print(f"  HUMAN SAFETY PROFILE v2 — GAS1 + BST2")
print(f"  Threshold: >{PRIMARY_THRESHOLD} | Patient-stratified | {len(donor_list)} donors")
print("="*70)
print(f"  {'Tissue':<38} {'Mean%':>7} {'±Std':>6}  {'Donors':>6}  Flag")
print("  " + "-"*65)
for _, row in report_df.iterrows():
    flag = "🚨 CRIT" if row["safety_flag"] == "CRITICAL" else \
           "⚠️  WARN" if row["safety_flag"] == "WARNING"  else "✓  SAFE"
    print(f"  {row['tissue']:<38} {row['coexpr_mean_pct']:>7.3f} {row['coexpr_std_pct']:>6.3f}  {row['n_donors']:>6}  {flag}")
print("="*70)

# ── PLOT: TRANSFER CURVES ─────────────────────────────────────────────────────
print("\nGenerating transfer curve plots...")

BG     = "#0e0e0f"
SURF   = "#161617"
PINK   = "#ff4d6d"
GREEN  = "#52d48a"
BLUE   = "#74b3ff"
YELLOW = "#ffd166"
MUTED  = "#888890"

# Pick tissues to highlight: critical organs + top co-expressors
highlight_tissues = list(report_df.head(12)["tissue"].values)

plt.style.use("dark_background")
fig = plt.figure(figsize=(20, 14), facecolor=BG)
gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)

ax1 = fig.add_subplot(gs[0, :])   # transfer curves — top tissues
ax2 = fig.add_subplot(gs[1, 0])   # bar chart at t=1.0
ax3 = fig.add_subplot(gs[1, 1])   # heatmap of all tissues x thresholds

fig.suptitle("GAS1 + BST2 Safety Profile — Patient-Stratified Transfer Curves\nAriadne Therapeutics",
             fontsize=13, color="white", fontweight="bold")

for ax in [ax1, ax2, ax3]:
    ax.set_facecolor(SURF)
    for spine in ax.spines.values():
        spine.set_edgecolor("#2a2a2d")
    ax.tick_params(colors=MUTED)

# ── PANEL 1: Transfer curves ──────────────────────────────────────────────────
ax1.set_facecolor(SURF)
cmap = plt.cm.get_cmap("tab20", len(highlight_tissues))

for i, tissue in enumerate(highlight_tissues):
    rates = [strat_results[tissue].get(t, {}).get("mean", 0.0) for t in THRESHOLDS]
    stds  = [strat_results[tissue].get(t, {}).get("std",  0.0) for t in THRESHOLDS]

    is_crit = any(o.lower() in tissue.lower() for o in CRITICAL_ORGANS)
    color   = PINK if is_crit else cmap(i)
    lw      = 2.5 if is_crit else 1.5
    alpha   = 1.0 if is_crit else 0.8

    ax1.plot(THRESHOLDS, rates, color=color, linewidth=lw, alpha=alpha,
             marker="o", markersize=4, label=tissue)
    ax1.fill_between(THRESHOLDS,
                     [r - s for r, s in zip(rates, stds)],
                     [r + s for r, s in zip(rates, stds)],
                     color=color, alpha=0.08)

ax1.axvline(PRIMARY_THRESHOLD, color=YELLOW, linestyle="--", linewidth=1.2,
            alpha=0.7, label=f"Primary threshold ({PRIMARY_THRESHOLD})")
ax1.axhline(5, color=PINK, linestyle=":", linewidth=1, alpha=0.5, label="5% safety line")
ax1.set_xlabel("Expression threshold (normalised counts per 10k)", color=MUTED, fontsize=10)
ax1.set_ylabel("Co-expression rate % (patient-stratified mean ± std)", color=MUTED, fontsize=10)
ax1.set_title("Transfer Curves — Co-expression rate vs threshold (top 12 tissues)", color="white", fontsize=11)
ax1.legend(fontsize=7.5, ncol=3, loc="upper right",
           facecolor=SURF, edgecolor="#2a2a2d", labelcolor="white")
ax1.grid(alpha=0.08)
ax1.set_xticks(THRESHOLDS)

# ── PANEL 2: Bar chart at primary threshold ───────────────────────────────────
plot_df = report_df.head(25).sort_values("coexpr_mean_pct", ascending=True)
bar_colors = [PINK if any(o.lower() in t.lower() for o in CRITICAL_ORGANS)
              else YELLOW if v > 5 else GREEN
              for t, v in zip(plot_df["tissue"], plot_df["coexpr_mean_pct"])]

bars = ax2.barh(range(len(plot_df)), plot_df["coexpr_mean_pct"],
                color=bar_colors, edgecolor="none", height=0.7)

# Error bars
ax2.errorbar(plot_df["coexpr_mean_pct"], range(len(plot_df)),
             xerr=plot_df["coexpr_std_pct"], fmt="none",
             color=MUTED, capsize=2, linewidth=1, alpha=0.6)

ax2.set_yticks(range(len(plot_df)))
ax2.set_yticklabels(plot_df["tissue"], fontsize=7.5)
ax2.set_xlabel("Co-expression rate % ± std", color=MUTED, fontsize=9)
ax2.set_title(f"Patient-Stratified Co-expression at threshold >{PRIMARY_THRESHOLD}\n(top 25 tissues)",
              color="white", fontsize=10)
ax2.axvline(5, color=PINK, linestyle="--", alpha=0.5, linewidth=1)
ax2.grid(alpha=0.08, axis="x")

# ── PANEL 3: Heatmap ──────────────────────────────────────────────────────────
# Top 20 tissues by max co-expression across all thresholds
top20 = sorted(tissue_list,
               key=lambda t: max(strat_results[t].get(th, {}).get("mean", 0) for th in THRESHOLDS),
               reverse=True)[:20]

heatmap_data = np.array([
    [strat_results[t].get(th, {}).get("mean", 0.0) for th in THRESHOLDS]
    for t in top20
])

im = ax3.imshow(heatmap_data, aspect="auto", cmap="RdYlGn_r",
                vmin=0, vmax=max(15, heatmap_data.max()))
ax3.set_xticks(range(len(THRESHOLDS)))
ax3.set_xticklabels([f">{t}" for t in THRESHOLDS], fontsize=8, color=MUTED)
ax3.set_yticks(range(len(top20)))
ax3.set_yticklabels(top20, fontsize=7.5, color="white")
ax3.set_xlabel("Threshold", color=MUTED, fontsize=9)
ax3.set_title("Co-expression Heatmap\n(top 20 tissues across all thresholds)", color="white", fontsize=10)

cbar = plt.colorbar(im, ax=ax3, fraction=0.03, pad=0.04)
cbar.set_label("Co-expression %", color=MUTED, fontsize=8)
cbar.ax.tick_params(colors=MUTED, labelsize=7)

# Annotate cells
for i in range(len(top20)):
    for j in range(len(THRESHOLDS)):
        val = heatmap_data[i, j]
        ax3.text(j, i, f"{val:.1f}", ha="center", va="center",
                 fontsize=6, color="white" if val > 8 else "#333")

plt.tight_layout()
out_path = os.path.join(OUT_DIR, "transfer_curves.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=BG)
plt.close()
print(f"Saved → {out_path}")
print(f"Saved → {OUT_DIR}/safety_report_t1.csv")
print(f"Saved → {OUT_DIR}/patient_stratified.csv")
print("\nDone. Send transfer_curves.png and safety_report_t1.csv to Haziq.")