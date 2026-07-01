"""
HuKSA — cluster-free kinase read-across predictor.

Replaces the HDBSCAN colony-label readout (shown to be unreliable on out-of-sample
compounds) with direct k-NN read-across over the 794 x 464 MEASURED bioactivity matrix:

  query SMILES
    -> Morgan/ECFP4 fingerprint
    -> Tanimoto similarity to all 794 reference compounds
    -> (1) TARGET    : Tanimoto-weighted mean measured pActivity over the k nearest
                       neighbours, per kinase (min support), ranked -> top-N targets
       (2) SELECTIVITY: Tanimoto-weighted read-across S-score over nearest neighbours -> tier
       (3) DOMAIN     : nearest-neighbour Tanimoto -> High/Moderate/Low/Outlier confidence

No UMAP, no HDBSCAN. The colony map (if wanted) is a separate *visualization* layer only.

CLI:
  python predictor.py "CCO..."            # JSON prediction for one SMILES
  python predictor.py --demo              # runs the worked examples (needs internet for PubChem)
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors as rd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")


def morgan_fp(smiles):
    m = Chem.MolFromSmiles(str(smiles))
    return rd.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) if m else None


class HuKSARx:
    """Cluster-free read-across predictor over the measured kinome matrix."""

    def __init__(self, data_dir=DATA):
        mat = pd.read_csv(os.path.join(data_dir, "compounds_x_kinases_pActivity.csv"))
        self.smiles_col = mat.columns[0]
        self.kin_cols = list(mat.columns[1:])                 # UniProt IDs
        self.M = mat[self.kin_cols].values.astype(float)      # 794 x 464 pActivity (NaN = unmeasured)
        self.smiles = mat[self.smiles_col].astype(str).tolist()
        self.fps = [morgan_fp(s) for s in self.smiles]
        meta = pd.read_csv(os.path.join(data_dir, "kinase_metadata_converted.csv"))
        gcol = "HGNCName" if "HGNCName" in meta.columns else "Name"
        ucol = [c for c in meta.columns if "uniprot" in c.lower()][0]
        u2g = {str(u): str(g) for u, g in zip(meta[ucol], meta[gcol])}
        self.col_gene = [u2g.get(u, u) for u in self.kin_cols]
        self.genes = sorted({g for g in self.col_gene if isinstance(g, str) and g})  # kinase panel
        # per-compound S-score = fraction of measured kinome engaged at pActivity > 6
        meas = ~np.isnan(self.M)
        hits = (np.nan_to_num(self.M, nan=-99) > 6.0) & meas
        with np.errstate(divide="ignore", invalid="ignore"):
            self.s_score = np.where(meas.sum(1) > 0, hits.sum(1) / meas.sum(1), np.nan)
        self.n = len(self.smiles)
        # --- precomputation for fast batch screening (BulkTanimoto needs non-None fps) ---
        self._vidx = [i for i, f in enumerate(self.fps) if f is not None]
        self._vfps = [self.fps[i] for i in self._vidx]
        self._vM = self.M[self._vidx]
        self._vs = self.s_score[self._vidx]
        self._gene_cols = {}
        for j, g in enumerate(self.col_gene):
            self._gene_cols.setdefault(g, []).append(j)

    @staticmethod
    def _confidence(nn):
        # heuristic applicability-domain bands on nearest-neighbour Tanimoto
        if nn >= 0.40: return "High"
        if nn >= 0.30: return "Moderate"
        if nn >= 0.20: return "Low"
        return "Outlier"

    @staticmethod
    def _tier(s):
        if s != s: return "Unknown"
        return "Highly Selective" if s <= 0.05 else ("Selective" if s <= 0.15 else "Promiscuous")

    def screen_row(self, smiles, k_target=5, k_sscore=15, min_support=2, target_gene=None):
        """Fast vectorized single-compound screen -> compact row dict (None if unparseable)."""
        q = morgan_fp(smiles)
        if q is None:
            return None
        sims = np.asarray(DataStructs.BulkTanimotoSimilarity(q, self._vfps))
        order = np.argsort(-sims)
        nn = float(sims[order[0]])
        valid_s = order[~np.isnan(self._vs[order])][:k_sscore]
        w = sims[valid_s]
        s = float((w * self._vs[valid_s]).sum() / w.sum()) if w.sum() > 0 else float("nan")
        tidx = order[:k_target]; tw = sims[tidx]; sub = self._vM[tidx]        # k x n_kinase
        mask = ~np.isnan(sub)
        cnt = mask.sum(0)
        wcol = (tw[:, None] * mask).sum(0)
        psum = (tw[:, None] * np.where(mask, sub, 0.0)).sum(0)
        with np.errstate(divide="ignore", invalid="ignore"):
            score = np.where((cnt >= min_support) & (wcol > 0), psum / wcol, -np.inf)
        top_j = int(np.argmax(score)) if np.isfinite(score).any() else -1
        row = {"SMILES": smiles,
               "top_target": self.col_gene[top_j] if top_j >= 0 else "—",
               "S_score": round(s, 3) if s == s else None,
               "tier": self._tier(s),
               "confidence": self._confidence(nn),
               "NN_Tanimoto": round(nn, 3)}
        if target_gene:
            vals = [score[j] for j in self._gene_cols.get(target_gene, []) if np.isfinite(score[j])]
            row[f"pAct[{target_gene}]"] = round(max(vals), 2) if vals else None
        return row

    def predict(self, smiles, k_target=5, k_sscore=15, top_n=10, min_support=2):
        qfp = morgan_fp(smiles)
        if qfp is None:
            return {"query_smiles": smiles, "valid": False, "error": "unparseable SMILES"}
        sims = np.array([DataStructs.TanimotoSimilarity(qfp, f) if f is not None else 0.0 for f in self.fps])
        order = np.argsort(-sims)
        nn = float(sims[order[0]])
        conf = self._confidence(nn)

        # selectivity: weighted read-across S-score over nearest neighbours with a defined S-score
        sidx = [i for i in order if not np.isnan(self.s_score[i])][:k_sscore]
        w = sims[sidx]; sval = self.s_score[np.array(sidx)]
        s = float((w * sval).sum() / w.sum()) if w.sum() > 0 else float("nan")

        # targets: weighted mean measured pActivity over k nearest neighbours, per kinase
        tidx = order[:k_target]; tw = sims[tidx]; sub = self.M[tidx]
        targets = []
        for j in range(self.M.shape[1]):
            col = sub[:, j]; ok = ~np.isnan(col)
            if ok.sum() >= min_support and tw[ok].sum() > 0:
                score = float((tw[ok] * col[ok]).sum() / tw[ok].sum())
                targets.append((self.col_gene[j], round(score, 2), int(ok.sum())))
        targets.sort(key=lambda x: -x[1])

        return {
            "query_smiles": smiles, "valid": True,
            "applicability_domain": {
                "nearest_neighbour_tanimoto": round(nn, 3),
                "confidence": conf,
                "note": "predictions below 'Low' (NN<0.30) are out-of-domain and weakly supported",
            },
            "selectivity": {"s_score": round(s, 3) if s == s else None, "tier": self._tier(s)},
            "predicted_targets": [
                {"rank": i + 1, "gene": g, "pred_pActivity": sc, "neighbours_measured": n}
                for i, (g, sc, n) in enumerate(targets[:top_n])
            ],
        }


# --------------------------------------------------------------------------- demo
# verified out-of-sample examples (SMILES from PubChem; none are in the 794 reference set)
EXAMPLES = [  # name, known target, SMILES
    ("nazartinib",   "EGFR",        "CC1=NC=CC(=C1)C(=O)NC2=NC3=C(N2[C@@H]4CCCCN(C4)C(=O)/C=C/CN(C)C)C(=CC=C3)Cl"),
    ("KU-0063794",   "MTOR",        "C[C@@H]1CN(C[C@@H](O1)C)C2=NC3=C(C=CC(=N3)C4=CC(=C(C=C4)OC)CO)C(=N2)N5CCOCC5"),
    ("tirabrutinib", "BTK",         "CC#CC(=O)N1CC[C@H](C1)N2C3=NC=NC(=C3N(C2=O)C4=CC=C(C=C4)OC5=CC=CC=C5)N"),
    ("sitravatinib", "MET/AXL/KDR", "COCCNCC1=CN=C(C=C1)C2=CC3=NC=CC(=C3S2)OC4=C(C=C(C=C4)NC(=O)C5(CC5)C(=O)NC6=CC=C(C=C6)F)F"),
    ("infigratinib", "FGFR1-4",     "CCN1CCN(CC1)C2=CC=C(C=C2)NC3=CC(=NC=N3)N(C)C(=O)NC4=C(C(=CC(=C4Cl)OC)OC)Cl"),
    ("selpercatinib","RET",         "CC(C)(COC1=CN2C(=C(C=N2)C#N)C(=C1)C3=CN=C(C=C3)N4CC5CC(C4)N5CC6=CN=C(C=C6)OC)O"),
]


def demo():
    rx = HuKSARx()
    print(f"loaded {rx.n} reference compounds x {len(rx.kin_cols)} kinases\n")
    for name, known, smi in EXAMPLES:
        r = rx.predict(smi)
        ad = r["applicability_domain"]; sel = r["selectivity"]
        top = ", ".join(f"{t['gene']}({t['pred_pActivity']})" for t in r["predicted_targets"][:5])
        hit = any(known.split("/")[0][:4] in t["gene"] or t["gene"] in known for t in r["predicted_targets"][:10])
        print(f"# {name}  (known: {known})")
        print(f"   domain      : NN-Tanimoto {ad['nearest_neighbour_tanimoto']}  -> confidence {ad['confidence']}")
        print(f"   selectivity : S={sel['s_score']}  -> {sel['tier']}")
        print(f"   top-5 targets: {top}")
        print(f"   known target in top-10? {hit}\n")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--demo":
        demo()
    elif len(sys.argv) > 1:
        print(json.dumps(HuKSARx().predict(sys.argv[1]), indent=2))
    else:
        print(__doc__)
