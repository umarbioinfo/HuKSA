"""
Leave-one-compound-out cross-validation of HuKSA-RX read-across over the 794x464 measured matrix.

For each reference compound i (self-excluded), predict its pActivity at every kinase it was
MEASURED at, from the Tanimoto-weighted mean of its k nearest OTHER compounds' measured activity
(min_support>=2). Score predicted-vs-measured with standard supervised metrics:

  Regression   : RMSE, MAE, R^2, Pearson, Spearman   (continuous pActivity)
  Classification (hit = pAct>6): pooled ROC-AUC, PR-AUC; per-compound AUC (rank kinases);
                 per-kinase AUC (rank compounds); recall@1/3/10 (target in top-N predicted)

CIs are compound-level bootstrap (resample compounds, not cells: cells within a compound are
correlated). No duplicate InChIKeys exist in the 794, so index self-exclusion is leakage-free.

Run from project root:  python manuscript/validation/tier1/loco/loco_cv.py
"""
import os, json, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors as rd
from sklearn.metrics import roc_auc_score, average_precision_score
from scipy.stats import pearsonr, spearmanr

OUT = "manuscript/validation/tier1/loco"; os.makedirs(OUT, exist_ok=True)
P   = "huksa-rx/app/data/compounds_x_kinases_pActivity.csv"   # the DEPLOYED matrix
HIT = 6.0; MIN_SUPPORT = 2; KS = [5, 9, 15]; KHEAD = 5; B = 1000
rng = np.random.default_rng(20260617)

mat = pd.read_csv(P); scol = mat.columns[0]; kin = list(mat.columns[1:])
M = mat[kin].values.astype(float)                       # 794 x 464, NaN = unmeasured
n, J = M.shape
def fp(s):
    m = Chem.MolFromSmiles(str(s)); return rd.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) if m else None
fps = [fp(s) for s in mat[scol]]
meas = ~np.isnan(M)
print(f"matrix {M.shape} | measured {meas.sum()} | hit prevalence {100*((M>HIT)&meas).sum()/meas.sum():.2f}%")

# ---- LOCO prediction: store cell-level (pred,true) per compound, for each k ----
def loco(k):
    per = []            # list of dicts per compound: idx, kin_idx, pred, true
    for i in range(n):
        if fps[i] is None: continue
        sims = np.array(DataStructs.BulkTanimotoSimilarity(fps[i], fps), dtype=float)
        sims[i] = -1.0                                   # SELF-EXCLUSION (no duplicates exist)
        order = np.argsort(-sims)[:k]; w = sims[order]
        sub = M[order]                                   # k x J  neighbour measured activity
        mk = ~np.isnan(sub); cnt = mk.sum(0)
        wsum = (w[:, None] * mk).sum(0)
        psum = (w[:, None] * np.where(mk, sub, 0.0)).sum(0)
        with np.errstate(divide="ignore", invalid="ignore"):
            pred = np.where((cnt >= MIN_SUPPORT) & (wsum > 0), psum / wsum, np.nan)
        ev = meas[i] & ~np.isnan(pred)                   # evaluate where i measured AND prediction exists
        kk = np.where(ev)[0]
        if kk.size:
            per.append(dict(idx=i, j=kk, pred=pred[kk], true=M[i, kk]))
    return per

# ---- metric helpers ----
def pooled(per):
    p = np.concatenate([d["pred"] for d in per]); t = np.concatenate([d["true"] for d in per])
    y = (t > HIT).astype(int)
    out = dict(n_cells=int(p.size), rmse=float(np.sqrt(np.mean((p-t)**2))), mae=float(np.mean(np.abs(p-t))),
               r2=float(1 - np.sum((p-t)**2)/np.sum((t-t.mean())**2)),
               pearson=float(pearsonr(p, t)[0]), spearman=float(spearmanr(p, t)[0]),
               roc_auc=float(roc_auc_score(y, p)) if 0 < y.mean() < 1 else float("nan"),
               pr_auc=float(average_precision_score(y, p)) if y.sum() > 0 else float("nan"),
               prevalence=float(y.mean()))
    return out

def per_compound_auc(per):
    aucs, recalls = [], {1: [], 3: [], 10: []}
    for d in per:
        y = (d["true"] > HIT).astype(int); s = d["pred"]
        if 0 < y.sum() < y.size:
            aucs.append(roc_auc_score(y, s))
        if y.sum() >= 1:                                  # recall@N: a true hit-kinase in top-N predicted
            ordr = np.argsort(-s)
            for N in recalls:
                recalls[N].append(int(y[ordr[:N]].sum() > 0))
    return np.array(aucs), {N: np.array(v) for N, v in recalls.items()}

def per_kinase_auc(per, min_comp=10, min_pos=3):
    by = {}
    for d in per:
        for j, pr, tr in zip(d["j"], d["pred"], d["true"]):
            by.setdefault(int(j), []).append((pr, tr))
    rows = []
    for j, vals in by.items():
        s = np.array([v[0] for v in vals]); y = (np.array([v[1] for v in vals]) > HIT).astype(int)
        if y.size >= min_comp and y.sum() >= min_pos and y.sum() < y.size:
            rows.append((j, kin[j], y.size, int(y.sum()), float(roc_auc_score(y, s))))
    return pd.DataFrame(rows, columns=["col", "uniprot", "n_comp", "n_pos", "auc"])

# ---- run headline k, store, and a k-robustness sweep ----
print("\n=== leave-one-compound-out CV ===")
sweep = []
per_head = None
for k in KS:
    per = loco(k)
    pm = pooled(per); pca, rec = per_compound_auc(per)
    cov = pm["n_cells"] / int(meas.sum())
    sweep.append(dict(k=k, coverage=round(cov, 3), rmse=round(pm["rmse"], 3), r2=round(pm["r2"], 3),
                      spearman=round(pm["spearman"], 3), pooled_auc=round(pm["roc_auc"], 3),
                      pr_auc=round(pm["pr_auc"], 3), per_cmpd_auc=round(float(np.mean(pca)), 3),
                      recall10=round(float(np.mean(rec[10])), 3)))
    print(f"k={k:2d} | cov {cov:.2f} | RMSE {pm['rmse']:.2f} R2 {pm['r2']:.2f} rho {pm['spearman']:.2f} "
          f"| pooledAUC {pm['roc_auc']:.3f} PR-AUC {pm['pr_auc']:.3f} | perCmpdAUC {np.mean(pca):.3f} "
          f"| recall@1/3/10 {np.mean(rec[1]):.2f}/{np.mean(rec[3]):.2f}/{np.mean(rec[10]):.2f}")
    if k == KHEAD:
        per_head = per; pm_head = pm; pca_head = pca; rec_head = rec
pd.DataFrame(sweep).to_csv(f"{OUT}/loco_k_sweep.csv", index=False)

# ---- per-kinase AUC at headline k ----
pk = per_kinase_auc(per_head); pk = pk.sort_values("auc", ascending=False)
pk.to_csv(f"{OUT}/loco_per_kinase_auc.csv", index=False)
print(f"\nper-kinase AUC (k={KHEAD}, kinases with >=10 compounds & >=3 hits, n={len(pk)}): "
      f"median {pk.auc.median():.3f} | mean {pk.auc.mean():.3f} | IQR [{pk.auc.quantile(.25):.3f},{pk.auc.quantile(.75):.3f}]")
print("  strongest:", ", ".join(f"{r.uniprot}({r.auc:.2f})" for _, r in pk.head(5).iterrows()))
print("  weakest  :", ", ".join(f"{r.uniprot}({r.auc:.2f})" for _, r in pk.tail(5).iterrows()))

# ---- compound-level bootstrap CIs (headline k) ----
def boot_ci(fn, B=B):
    idx = np.arange(len(per_head)); stats = []
    for _ in range(B):
        samp = rng.choice(idx, size=idx.size, replace=True)
        stats.append(fn([per_head[s] for s in samp]))
    lo, hi = np.nanpercentile(stats, [2.5, 97.5]); return float(lo), float(hi)
ci_rmse = boot_ci(lambda pr: pooled(pr)["rmse"])
ci_auc  = boot_ci(lambda pr: pooled(pr)["roc_auc"])
ci_pr   = boot_ci(lambda pr: pooled(pr)["pr_auc"])
ci_pca  = boot_ci(lambda pr: float(np.mean(per_compound_auc(pr)[0])))
ci_r10  = boot_ci(lambda pr: float(np.mean(per_compound_auc(pr)[1][10])))

# ---- cell-level predictions (headline k) for the record ----
recs = []
for d in per_head:
    for j, pr, tr in zip(d["j"], d["pred"], d["true"]):
        recs.append((mat[scol].iloc[d["idx"]], kin[j], round(float(pr), 3), round(float(tr), 3), int(tr > HIT)))
pd.DataFrame(recs, columns=["smiles", "uniprot", "pred_pAct", "meas_pAct", "is_hit"]).to_csv(
    f"{OUT}/loco_cell_predictions_k{KHEAD}.csv", index=False)

summary = dict(k=KHEAD, n_compounds=len(per_head), coverage=round(pm_head["n_cells"]/int(meas.sum()), 3),
   n_cells=pm_head["n_cells"], prevalence=round(pm_head["prevalence"], 4),
   rmse=[round(pm_head["rmse"],3), [round(ci_rmse[0],3), round(ci_rmse[1],3)]],
   r2=round(pm_head["r2"],3), spearman=round(pm_head["spearman"],3),
   pooled_roc_auc=[round(pm_head["roc_auc"],3), [round(ci_auc[0],3), round(ci_auc[1],3)]],
   pr_auc=[round(pm_head["pr_auc"],3), [round(ci_pr[0],3), round(ci_pr[1],3)]],
   per_compound_auc=[round(float(np.mean(pca_head)),3), [round(ci_pca[0],3), round(ci_pca[1],3)]],
   recall_at_1=round(float(np.mean(rec_head[1])),3), recall_at_3=round(float(np.mean(rec_head[3])),3),
   recall_at_10=[round(float(np.mean(rec_head[10])),3), [round(ci_r10[0],3), round(ci_r10[1],3)]],
   per_kinase_auc_median=round(float(pk.auc.median()),3), per_kinase_auc_n=int(len(pk)))
json.dump(summary, open(f"{OUT}/loco_summary.json", "w"), indent=2)

print("\n================  HEADLINE (k=5)  ================")
print(f"compounds evaluated      : {summary['n_compounds']}/{n}")
print(f"cells predicted          : {summary['n_cells']}  (coverage {summary['coverage']*100:.0f}% of measured)")
print(f"hit prevalence           : {summary['prevalence']*100:.1f}%")
print(f"RMSE (pActivity)         : {summary['rmse'][0]}  95%CI {summary['rmse'][1]}")
print(f"R^2 / Spearman           : {summary['r2']} / {summary['spearman']}")
print(f"pooled ROC-AUC           : {summary['pooled_roc_auc'][0]}  95%CI {summary['pooled_roc_auc'][1]}")
print(f"PR-AUC (base {summary['prevalence']*100:.1f}%)   : {summary['pr_auc'][0]}  95%CI {summary['pr_auc'][1]}")
print(f"per-compound AUC (rank kinases) : {summary['per_compound_auc'][0]}  95%CI {summary['per_compound_auc'][1]}")
print(f"per-kinase  AUC (rank compounds) median : {summary['per_kinase_auc_median']}  (n={summary['per_kinase_auc_n']} kinases)")
print(f"recall@1/3/10 (target in top-N) : {summary['recall_at_1']} / {summary['recall_at_3']} / {summary['recall_at_10'][0]}  (top10 95%CI {summary['recall_at_10'][1]})")
print("saved ->", OUT)
