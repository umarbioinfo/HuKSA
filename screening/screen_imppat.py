# -*- coding: utf-8 -*-
"""
Screen the IMPPAT phytochemical library through HuKSA.

Reads every *_3D.sdf in IMPPAT LIBRARY/LIbrary/, fingerprints it, and runs the
HuKSA cluster-free k-NN screen (screen_row) to get, per compound:
  IMPHY id, SMILES, predicted top kinase target, S-score, selectivity tier,
  applicability-domain confidence band, and nearest-neighbour Tanimoto.

The headline result is the APPLICABILITY-DOMAIN / QUARANTINE breakdown: how much of
a natural-product library falls OUTSIDE the chemical space HuKSA was built on, and
is therefore correctly flagged as out-of-domain rather than given false confidence.

Run:  python huksa-rx/IMPPAT/screen_imppat.py
Out:  huksa-rx/IMPPAT/imppat_huksa_screen.csv
"""
import os, sys, glob, time
HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.join(HERE, "..", "app")
sys.path.insert(0, APP)
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
from predictor import HuKSARx

LIB = os.path.join(HERE, "IMPPAT LIBRARY", "LIbrary")

def main():
    rx = HuKSARx()
    files = sorted(glob.glob(os.path.join(LIB, "*.sdf")) + glob.glob(os.path.join(LIB, "*.mol")))
    n_sdf = sum(f.lower().endswith(".sdf") for f in files)
    n_mol = sum(f.lower().endswith(".mol") for f in files)
    print(f"reference matrix: {rx.n} compounds x {len(rx.kin_cols)} kinases")
    print(f"IMPPAT structure files found: {len(files)}  ({n_sdf} .sdf + {n_mol} .mol)")

    rows, bad_parse, t0 = [], 0, time.time()
    for n, f in enumerate(files, 1):
        imphy = os.path.basename(f).replace("_3D.sdf", "").replace("_3D.mol", "")
        try:
            if f.lower().endswith(".mol"):
                m = Chem.MolFromMolFile(f, removeHs=False)
            else:
                suppl = Chem.SDMolSupplier(f, removeHs=False)
                m = suppl[0] if len(suppl) else None
        except Exception:
            m = None
        if m is None:
            bad_parse += 1
            continue
        smi = Chem.MolToSmiles(m)
        r = rx.screen_row(smi)
        if r is None:
            bad_parse += 1
            continue
        r = {"IMPHY_ID": imphy, **r}
        rows.append(r)
        if n % 1000 == 0:
            print(f"  ...{n}/{len(files)}  ({time.time()-t0:.0f}s)")

    df = pd.DataFrame(rows)
    out = os.path.join(HERE, "imppat_huksa_screen.csv")
    df.to_csv(out, index=False)
    print(f"\nwrote {out}")
    print(f"screened ok: {len(df)}   unparseable: {bad_parse}   time: {time.time()-t0:.0f}s\n")

    # ---- the quarantine / applicability-domain story ----
    order = ["High", "Moderate", "Low", "Outlier"]
    cb = df["confidence"].value_counts().reindex(order).fillna(0).astype(int)
    pct = (cb / len(df) * 100).round(1)
    print("APPLICABILITY-DOMAIN (confidence) breakdown:")
    for b in order:
        print(f"   {b:9s}: {cb[b]:6d}  ({pct[b]:5.1f}%)")
    out_of_domain = int(cb["Low"] + cb["Outlier"])
    print(f"   --> OUT-OF-DOMAIN (Low+Outlier, NN<0.30): {out_of_domain} ({out_of_domain/len(df)*100:.1f}%)")
    in_domain = int(cb["High"] + cb["Moderate"])
    print(f"   --> IN-DOMAIN     (High+Moderate, NN>=0.30): {in_domain} ({in_domain/len(df)*100:.1f}%)\n")

    print("Selectivity tier among IN-DOMAIN compounds (NN>=0.30):")
    ind = df[df["confidence"].isin(["High", "Moderate"])]
    if len(ind):
        print(ind["tier"].value_counts().to_string())
        print("\nTop predicted kinase targets among IN-DOMAIN compounds:")
        print(ind["top_target"].value_counts().head(12).to_string())
    print(f"\nNN-Tanimoto distribution: median {df['NN_Tanimoto'].median():.3f}  "
          f"max {df['NN_Tanimoto'].max():.3f}  "
          f"95th pct {df['NN_Tanimoto'].quantile(0.95):.3f}")

if __name__ == "__main__":
    main()
