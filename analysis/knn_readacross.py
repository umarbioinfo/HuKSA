"""
k-NN-over-measured-profiles target prediction, as a head-to-head
against HuKSA's colony-inheritance. For each out-of-sample query, predict targets directly from
the Tanimoto-weighted MEASURED kinase activity of its nearest neighbours in the 794-set.
Run from project root:  python manuscript/validation/tier1/knn/knn_readacross.py
"""
import os, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")
from rdkit import Chem, DataStructs
from rdkit.Chem import rdMolDescriptors as rd

OUT = "manuscript/validation/tier1/knn"; os.makedirs(OUT, exist_ok=True)
D = "huksa-atlas/app/saved_atlas_tool/"

# ---- combined OOS query set (from the two frozen runs) + acceptable target genes (family-aware) ----
TGT = {
 "mobocertinib":{"EGFR"},"olmutinib":{"EGFR"},"nazartinib":{"EGFR"},"brigatinib":{"ALK"},
 "lorlatinib":{"ALK","ROS1"},"radotinib":{"ABL1"},
 "nintedanib":{"KDR","FGFR1","FGFR2","FGFR3","PDGFRA","PDGFRB","FLT1","FLT4"},
 "infigratinib":{"FGFR1","FGFR2","FGFR3","FGFR4"},"erdafitinib":{"FGFR1","FGFR2","FGFR3","FGFR4"},
 "pemigatinib":{"FGFR1","FGFR2","FGFR3"},"upadacitinib":{"JAK1"},"acalabrutinib":{"BTK"},
 "tirabrutinib":{"BTK"},"zanubrutinib":{"BTK"},"pirtobrutinib":{"BTK"},"asciminib":{"ABL1"},
 "selpercatinib":{"RET"},"repotrectinib":{"ROS1","NTRK1","NTRK2","NTRK3","ALK"},
 "SP600125":{"MAPK8","MAPK9","MAPK10"},"U0126":{"MAP2K1","MAP2K2"},"PD173074":{"FGFR1","FGFR2","FGFR3"},
 "SU11274":{"MET"},"GSK429286A":{"ROCK1","ROCK2"},"H-89":{"PRKACA","PRKACB"},"KU-0063794":{"MTOR"},
 "Torin 2":{"MTOR"},"LY294002":{"PIK3CA","PIK3CB","PIK3CD","PIK3CG","MTOR"},"SCH772984":{"MAPK1","MAPK3"},
 "dorsomorphin":{"PRKAA1","PRKAA2"},"sitravatinib":{"MET","AXL","KDR","FLT1","FLT4"},
}
ALIAS = {"FMS":"CSF1R","PKACA":"PRKACA","PKACB":"PRKACB","AMPKA1":"PRKAA1","AMPKA2":"PRKAA2","CHK1":"CHEK1",
 "CHK2":"CHEK2","FRAP":"MTOR","FRAP1":"MTOR","VEGFR2":"KDR","VEGFR1":"FLT1","VEGFR3":"FLT4","MEK1":"MAP2K1",
 "MEK2":"MAP2K2","ERK1":"MAPK3","ERK2":"MAPK1","JNK1":"MAPK8","JNK2":"MAPK9","JNK3":"MAPK10","P38A":"MAPK14",
 "P38":"MAPK14","ABL":"ABL1","AURA":"AURKA","AURB":"AURKB","PDGFRA":"PDGFRA","AXL":"AXL"}
def nrm(g):
    g=str(g).upper()
    for s in ("_B","_A"):
        if g.endswith(s): g=g[:-2]
    return ALIAS.get(g,g)

def fp(s):
    m=Chem.MolFromSmiles(str(s)); return rd.GetMorganFingerprintAsBitVect(m,2,nBits=2048) if m else None

# query SMILES from the two frozen prediction files
q1=pd.read_csv("manuscript/validation/tier1/huksa_predictions_FROZEN.csv")
q2=pd.read_csv("manuscript/validation/tier1/toolcompounds/huksa_predictions_FROZEN.csv")
queries=pd.concat([q1[["name","smiles"]],q2[["name","smiles"]]]).drop_duplicates("name").reset_index(drop=True)
queries=queries[queries.name.isin(TGT)].reset_index(drop=True)

# measured matrix + uniprot->gene
mat=pd.read_csv(D+"compounds_x_kinases_pActivity.csv")
smiles_col=mat.columns[0]
kin_cols=list(mat.columns[1:])                       # UniProt IDs
M=mat[kin_cols].values.astype(float)                 # 794 x 464 pActivity (NaN = unmeasured)
atlas_fps=[fp(s) for s in mat[smiles_col]]
meta=pd.read_csv("saved_atlas_tool/kinase_metadata_converted.csv")
gcol="HGNCName" if "HGNCName" in meta.columns else "Name"   # HGNC = standard symbol
ucol=[c for c in meta.columns if "uniprot" in c.lower()][0]
uni2gene={str(u):nrm(g) for u,g in zip(meta[ucol],meta[gcol])}
col_gene=[uni2gene.get(u,nrm(u)) for u in kin_cols]   # normalized gene per matrix column
print(f"matrix {M.shape} | queries {len(queries)} | kinase cols mapped to genes: {sum(1 for g in col_gene if g)}")

def knn_predict(qsmi, k, min_meas=2):
    qfp=fp(qsmi)
    sims=np.array([DataStructs.TanimotoSimilarity(qfp,f) if f is not None else 0.0 for f in atlas_fps])
    idx=np.argsort(-sims)[:k]; w=sims[idx]
    sub=M[idx]                                        # k x 464
    score=np.full(M.shape[1], -np.inf); meas_n=np.zeros(M.shape[1],int)
    for j in range(M.shape[1]):
        col=sub[:,j]; ok=~np.isnan(col)
        meas_n[j]=ok.sum()
        if ok.sum()>=min_meas and w[ok].sum()>0:
            score[j]=(w[ok]*col[ok]).sum()/w[ok].sum()
    order=np.argsort(-score)
    ranked=[(col_gene[j], round(float(score[j]),2)) for j in order if np.isfinite(score[j])]
    return float(sims[idx].max()), ranked, meas_n

def hit_at(ranked, acc, n):
    accn={nrm(a) for a in acc}
    return any(g in accn for g,_ in ranked[:n])

# sweep k
print("\n=== k-NN read-across accuracy (known target in predicted top-N) ===")
rows=[]
for k in (5,9,15,30):
    res={"k":k}
    for n in (1,3,5,10):
        hits=0; total=0
        for _,r in queries.iterrows():
            _,ranked,_=knn_predict(r.smiles,k)
            total+=1; hits+=hit_at(ranked,TGT[r["name"]],n)
        res[f"top{n}"]=round(100*hits/total,1)
    rows.append(res); print(res)
pd.DataFrame(rows).to_csv(f"{OUT}/knn_accuracy_sweep.csv",index=False)

# detailed per-compound at k=9
K=9
det=[]
for _,r in queries.iterrows():
    nn,ranked,meas_n=knn_predict(r.smiles,K)
    accn={nrm(a) for a in TGT[r["name"]]}
    # is the true target even measured among the k neighbours?
    measurable=any((col_gene[j] in accn) and (meas_n[j]>=2) for j in range(M.shape[1]))
    top5=[g for g,_ in ranked[:5]]
    det.append(dict(name=r["name"], nn_tanimoto=round(nn,2), known=";".join(sorted(TGT[r["name"]])),
                    knn_top5="; ".join(top5), hit_top1=hit_at(ranked,TGT[r["name"]],1),
                    hit_top3=hit_at(ranked,TGT[r["name"]],3), hit_top10=hit_at(ranked,TGT[r["name"]],10),
                    target_measurable=measurable))
det=pd.DataFrame(det); det.to_csv(f"{OUT}/knn_per_compound_k9.csv",index=False)
pd.set_option("display.width",230); pd.set_option("display.max_colwidth",40)
print(f"\n=== per-compound (k={K}) ===")
print(det[["name","nn_tanimoto","known","knn_top5","hit_top3","hit_top10","target_measurable"]].to_string(index=False))

print("\n=== summary (k=9) ===")
print(f"top-1 {det.hit_top1.mean()*100:.0f}% | top-3 {det.hit_top3.mean()*100:.0f}% | top-10 {det.hit_top10.mean()*100:.0f}%")
print(f"target measurable among neighbours: {det.target_measurable.mean()*100:.0f}% ({det.target_measurable.sum()}/{len(det)})")
meas=det[det.target_measurable]
print(f"among compounds whose target IS measurable (n={len(meas)}): top-3 {meas.hit_top3.mean()*100:.0f}% | top-10 {meas.hit_top10.mean()*100:.0f}%")
print(f"\n(for comparison, colony-inheritance gave ~3% exact / ~12% family across the same compounds)")
