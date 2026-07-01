#!/usr/bin/env python3
"""
reproduce_selectivity_analysis.py
=================================================================
Reproducibility script for the HuKSA selectivity-metric analysis.

Regenerates EVERY quantitative result reported in the "Selectivity metric
selection" subsection, in labelled sections so each printed value maps to a
manuscript claim. Self-contained: depends only on pandas/numpy/scipy/matplotlib.

Run:
    python reproduce_selectivity_analysis.py

Inputs (paths relative to this script's location):
    saved_atlas_tool/compounds_x_kinases_pActivity.csv   794 x 464 pActivity matrix
    saved_atlas_tool/atlas_background_data.csv            SMILES -> Cluster_ID
    all_bioactivity_with_censoring.csv                    long-format raw (assay type, pActivity)

Outputs (written next to this script):
    repro_per_compound_metrics.csv     all per-compound metrics
    repro_cluster_verdicts.csv         per-cluster mean S-score + verdict
    repro_kd_validation.csv            Kd-subset validation table
    figure_selectivity_metric.png/.pdf 4-panel metric-selection figure

Expected key numbers (this dataset):
    r(Gini_padded, coverage)            = -0.997
    Gini measured-only median           ~ 0.03-0.05  (degenerate)
    entropy(log) median                 = 0.999       (degenerate)
    S-score(tau=6) coverage r           = -0.12
    S6 median / range                   = 0.095 / 0-0.95
    threshold stability S6 vs S7        = 0.85
    Kd validation set                   = 272 compounds (>=100 Kd kinases)
    S(pooled) vs S(Kd)                  = 0.97
    S(Kd) vs Kd-entropy (>=6/>=7/>=8)   = 0.37 / 0.61 / 0.81
=================================================================
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use("Agg")           # headless-safe
import matplotlib.pyplot as plt

HERE       = os.path.dirname(os.path.abspath(__file__))
OUT_DIR    = os.path.join(HERE, "repro_outputs")
os.makedirs(OUT_DIR, exist_ok=True)
# Reference data ships with this repo under app/data/. The large long-format raw file
# (all_bioactivity_with_censoring.csv, ~99 MB) is NOT bundled — drop it in analysis/data/
# to enable the optional Kd-validation section (Section 6); the script skips it if absent.
MATRIX_CSV = os.path.join(HERE, "..", "app", "data", "compounds_x_kinases_pActivity.csv")
ATLAS_CSV  = os.path.join(HERE, "..", "app", "data", "atlas_background_data.csv")
LONG_CSV   = os.path.join(HERE, "data", "all_bioactivity_with_censoring.csv")

THR      = 6.0      # hit threshold (pActivity > 6.0  == potency < 1 uM)
MIN_MEAS = 100      # panel-size floor (matches upstream pipeline filter)
MIN_KD   = 100      # min Kd-measured kinases for the validation subset


# ----------------------------------------------------------------------
# metric definitions
# ----------------------------------------------------------------------
def gini(values):
    """Standard sorted-values Gini on a non-negative 1-D array."""
    a = np.sort(np.asarray(values, dtype=float))
    n = a.size
    s = a.sum()
    if n == 0 or s == 0:
        return 0.0
    idx = np.arange(1, n + 1)
    return float((np.sum((2 * idx - n - 1) * a)) / (n * s))


def gini_padded(row):
    """Gini over the FULL panel, unmeasured -> 0 (the original, biased version)."""
    a = np.nan_to_num(row.values.astype(float), nan=0.0)
    return gini(np.clip(a, 0.0, None))


def gini_measured(row):
    """Gini over measured kinases only (clip negatives)."""
    m = row.dropna()
    if len(m) < MIN_MEAS:
        return np.nan
    return gini(np.clip(m.values.astype(float), 0.0, None))


def entropy_norm(values):
    """Normalized Shannon entropy of a non-negative profile."""
    v = np.clip(np.asarray(values, dtype=float), 0.0, None)
    s = v.sum()
    if s <= 0 or v.size < 2:
        return np.nan
    p = v[v > 0] / s
    return float(-(p * np.log(p)).sum() / np.log(v.size))


def entropy_log(row):
    """Entropy on log proportions (p_i ~ pActivity_i) -- the degenerate recipe."""
    m = row.dropna()
    if len(m) < MIN_MEAS:
        return np.nan
    return entropy_norm(m.values)


def entropy_lin(row):
    """Entropy on linear proportions (p_i ~ 10^pActivity_i)."""
    m = row.dropna()
    if len(m) < MIN_MEAS:
        return np.nan
    return entropy_norm(np.power(10.0, m.values.astype(float)))


def gini_linear(row):
    m = row.dropna()
    if len(m) < MIN_MEAS:
        return np.nan
    return gini(np.power(10.0, m.values.astype(float)))


def s_score(row, thr):
    m = row.dropna()
    if len(m) < MIN_MEAS:
        return np.nan
    return float((m > thr).sum()) / len(m)


def sel_window(row):
    m = row.dropna()
    if len(m) < MIN_MEAS:
        return np.nan
    return float(m.max() - m.median())


def banner(title):
    print("\n" + "=" * 70 + f"\n{title}\n" + "=" * 70)


# ----------------------------------------------------------------------
def main():
    df = pd.read_csv(MATRIX_CSV, index_col=0)
    coverage = df.notna().sum(axis=1)
    print(f"matrix: {df.shape[0]} compounds x {df.shape[1]} kinases  "
          f"({100*df.isna().sum().sum()/df.size:.1f}% missing)")

    m = pd.DataFrame(index=df.index)
    m["Measured_Count"] = coverage
    m["Gini_padded"]    = df.apply(gini_padded,   axis=1)
    m["Gini_measured"]  = df.apply(gini_measured, axis=1)
    m["Entropy_log"]    = df.apply(entropy_log,   axis=1)
    m["Entropy_lin"]    = df.apply(entropy_lin,   axis=1)
    m["Gini_linear"]    = df.apply(gini_linear,   axis=1)
    for t in (5.0, 6.0, 7.0):
        m[f"S{int(t)}"] = df.apply(lambda r, t=t: s_score(r, t), axis=1)
    m["S_Score"]    = m["S6"]
    m["Sel_Window"] = df.apply(sel_window, axis=1)
    v = m.dropna(subset=["S_Score"])

    # --- Section 1: Gini coverage confound ---
    banner("1. Gini is coverage-confounded")
    print(f"  r(Gini_padded, coverage)   = {pearsonr(coverage, m.Gini_padded)[0]:+.4f}   [report: -0.997]")
    print(f"  r(Gini_measured, coverage) = {pearsonr(v.Measured_Count, v.Gini_measured)[0]:+.4f}")

    # --- Section 2: collapse of measured-only Gini and log entropy ---
    banner("2. Degeneracy on log scale (measured-only)")
    print(f"  Gini_measured median = {v.Gini_measured.median():.4f}   [report: ~0.03-0.05]")
    print(f"  Entropy_log  median  = {v.Entropy_log.median():.4f}   IQR "
          f"[{v.Entropy_log.quantile(.25):.4f}, {v.Entropy_log.quantile(.75):.4f}]   [report: 0.999]")
    # toy example: max-selective vs fully-promiscuous on log proportions
    sel = np.array([9.0] + [5.0] * 49)
    pro = np.array([7.0] * 50)
    print(f"  toy log-entropy  selective={entropy_norm(sel):.4f}  promiscuous={entropy_norm(pro):.4f}  "
          f"[report: 0.9987 vs 1.000]")
    print(f"  toy lin-entropy  selective={entropy_norm(10**sel):.4f}  promiscuous={entropy_norm(10**pro):.4f}")

    # --- Section 3: linear restores discrimination; coverage bias stays mild ---
    #     (the obstacle to linearization is pooled-assay inconsistency, not coverage)
    banner("3. Linearization: discrimination returns; coverage bias stays mild")
    print(f"  Entropy_lin median = {v.Entropy_lin.median():.3f}  (vs 0.999 log)")
    print(f"  Gini_linear median = {v.Gini_linear.median():.3f}  (vs ~0.05 log)")
    print(f"  r(Gini_linear, coverage) = {pearsonr(v.Measured_Count, v.Gini_linear)[0]:+.3f}   [report: +0.16]")

    # --- Section 4: S-score properties ---
    banner("4. S-score: coverage-clean, discriminating, threshold-stable")
    print(f"  coverage r(S6) = {pearsonr(v.Measured_Count, v.S6)[0]:+.3f}   [report: -0.12]")
    print(f"  median S5/S6/S7 = {v.S5.median():.3f} / {v.S6.median():.3f} / {v.S7.median():.3f}")
    print(f"  S6 range = {v.S6.min():.3f} - {v.S6.max():.3f}")
    print(f"  rank stability  S6 vs S7 = {spearmanr(v.S6, v.S7).correlation:.3f}   "
          f"S6 vs S5 = {spearmanr(v.S6, v.S5).correlation:.3f}")

    # --- Section 5: cluster verdicts + face validity ---
    banner("5. Cluster verdicts and face validity")
    atlas = pd.read_csv(ATLAS_CSV)
    atlas = atlas.drop(columns=[c for c in ["S_Score"] if c in atlas.columns])  # avoid merge collision
    atlas = atlas.merge(v[["S_Score"]], left_on="SMILES", right_index=True, how="left")
    cl = atlas[atlas.Cluster_ID != -1].groupby("Cluster_ID")["S_Score"].mean().sort_values()
    # Fixed interpretable cutoffs: S = fraction of kinome engaged.
    # Highly Selective <=5%, Selective <=15%, Promiscuous >15% (coincide with dataset tertiles).
    t1, t2 = 0.05, 0.15
    verdict = cl.apply(lambda x: "Highly Selective" if x <= t1 else ("Selective" if x <= t2 else "Promiscuous"))
    print(f"  most selective cluster:  {cl.index[0]}  mean S = {cl.iloc[0]:.4f}  ({verdict.iloc[0]})")
    print(f"  most promiscuous cluster: {cl.index[-1]}  mean S = {cl.iloc[-1]:.4f}  ({verdict.iloc[-1]})")
    cl_out = pd.DataFrame({"Mean_S_Score": cl.round(4), "Verdict": verdict})

    # --- Section 6: Kd-subset validation (optional; needs the large raw long-format file) ---
    banner("6. Validation against gold-standard linear Kd entropy")
    if os.path.exists(LONG_CSV):
        kdf = kd_validation(v["S6"])
    else:
        print("  [skipped] analysis/data/all_bioactivity_with_censoring.csv not found -")
        print("            Section 6 (Kd validation) and figure panel C are omitted.")
        kdf = None

    # --- write outputs ---
    m.to_csv(os.path.join(OUT_DIR, "repro_per_compound_metrics.csv"))
    cl_out.to_csv(os.path.join(OUT_DIR, "repro_cluster_verdicts.csv"))
    if kdf is not None:
        kdf.to_csv(os.path.join(OUT_DIR, "repro_kd_validation.csv"), index=False)
    make_figure(coverage, m, v, cl, kdf)
    _kd = "repro_kd_validation.csv, " if kdf is not None else ""
    print(f"\nwrote: repro_per_compound_metrics.csv, repro_cluster_verdicts.csv, "
          f"{_kd}figure_selectivity_metric.png/.pdf")


def kd_validation(pooled_s6):
    cols = ["Standardized_SMILES", "UniProt ID", "Activity_Type_Std", "pActivity"]
    long = pd.read_csv(LONG_CSV, usecols=cols)
    kd = long[(long.Activity_Type_Std == "Kd") &
              (long.Standardized_SMILES.isin(pooled_s6.index))]
    pair = kd.groupby(["Standardized_SMILES", "UniProt ID"])["pActivity"].median().reset_index()
    rows = []
    for smiles, g in pair.groupby("Standardized_SMILES"):
        x = g.pActivity.values
        if len(x) < MIN_KD:
            continue
        lin = np.power(10.0, x)
        rows.append((smiles, len(x), float((x > THR).mean()),
                     float(-(lin/lin.sum() * np.log(lin/lin.sum())).sum() / np.log(len(x))),
                     float(x.max()), float(pooled_s6[smiles])))
    kdf = pd.DataFrame(rows, columns=["SMILES", "n_Kd", "S_Kd", "Kd_Entropy", "max_pKd", "S_pooled"])
    print(f"  validation set: {len(kdf)} compounds (mean panel {kdf.n_Kd.mean():.0f})")
    print(f"  S(pooled) vs S(Kd) = {spearmanr(kdf.S_pooled, kdf.S_Kd).correlation:+.3f}   [report: 0.97]")
    for f in (0, 6, 7, 8):
        s = kdf[kdf.max_pKd >= f]
        print(f"  S(Kd) vs Kd-entropy  max pKd>={f}: n={len(s):3d}  "
              f"rho={spearmanr(s.S_Kd, s.Kd_Entropy).correlation:+.3f}")
    return kdf


def make_figure(coverage, m, v, cl, kdf):
    from matplotlib.patches import Patch
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
                         "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9})
    fig, ax = plt.subplots(2, 2, figsize=(11, 9), layout="constrained")

    def tag(a, letter):                       # bold standalone panel letter
        a.text(-0.13, 1.04, letter, transform=a.transAxes, fontsize=14, fontweight="bold", va="top")

    # A: Gini confound
    ax[0,0].scatter(coverage, m.Gini_padded, s=8, alpha=.4)
    ax[0,0].set(xlabel="kinases measured", ylabel="Gini (padded)",
                title=f"Gini vs coverage (r = {pearsonr(coverage, m.Gini_padded)[0]:.3f})"); tag(ax[0,0], "A")
    # B: collapse
    ax[0,1].hist(v.Gini_measured.dropna(), bins=40, alpha=.6, label="Gini (log, measured)")
    ax[0,1].hist(v.Entropy_log.dropna(), bins=40, alpha=.6, label="Entropy (log)")
    ax[0,1].set(xlabel="metric value (0–1)", ylabel="compounds", title="Degeneracy on log scale")
    ax[0,1].legend(loc="upper right", framealpha=.9); tag(ax[0,1], "B")
    # C: validation (skipped if the raw Kd data was not bundled)
    if kdf is not None:
        sc = ax[1,0].scatter(kdf.S_Kd, kdf.Kd_Entropy, c=kdf.max_pKd, s=14, cmap="viridis")
        plt.colorbar(sc, ax=ax[1,0], label="max pKd")
        ax[1,0].set(xlabel="S-score (Kd)", ylabel="Kd entropy (linear)",
                    title="S-score vs gold-standard Kd entropy")
    else:
        ax[1,0].text(0.5, 0.5, "Kd validation skipped\n(raw file not bundled)",
                     ha="center", va="center", transform=ax[1,0].transAxes, color="grey")
        ax[1,0].set(xticks=[], yticks=[], title="S-score vs gold-standard Kd entropy")
    tag(ax[1,0], "C")
    # D: face validity — fixed interpretable cutoffs (Highly Selective <=0.05, Selective <=0.15)
    colors = ["#2c7fb8" if x <= 0.05 else ("#fec44f" if x <= 0.15 else "#d95f0e") for x in cl]
    ax[1,1].bar(range(len(cl)), cl.values, color=colors)
    ax[1,1].axhline(0.05, ls=":", c="grey", lw=.8); ax[1,1].axhline(0.15, ls=":", c="grey", lw=.8)
    ax[1,1].set(xlabel="cluster (sorted)", ylabel="mean S-score", title="Per-cluster selectivity")
    ax[1,1].legend(handles=[Patch(color="#2c7fb8", label="Highly Selective (≤0.05)"),
                            Patch(color="#fec44f", label="Selective (≤0.15)"),
                            Patch(color="#d95f0e", label="Promiscuous (>0.15)")],
                   loc="upper left", framealpha=.9); tag(ax[1,1], "D")
    fig.savefig(os.path.join(OUT_DIR, "figure_selectivity_metric.png"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(OUT_DIR, "figure_selectivity_metric.pdf"), bbox_inches="tight")


if __name__ == "__main__":
    main()
