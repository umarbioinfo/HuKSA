# HuKSA — Human Kinome Selectivity Atlas

A transparent, ligand-based tool that predicts a small molecule's likely **kinase targets**, its
**global selectivity**, and an explicit **confidence flag** — directly from a SMILES string, with no
model training. HuKSA is built around a simple principle: **know when *not* to trust the prediction.**
Every prediction exposes the measured reference compounds behind it and flags when a query falls
outside the chemistry the method can speak to.

**Live web app:** https://huksalab-tool.hf.space

---

## Method in brief

```
SMILES ──> 2048-bit Morgan/ECFP4 fingerprint ──> Tanimoto similarity to 794 reference inhibitors
       ──> Tanimoto-weighted k-nearest-neighbour (k = 5) read-across over the measured activity matrix
       ──> per-kinase predicted potency  +  ranked targets
       ──> S-score (fraction of the measured kinome engaged at pActivity > 6)  +  selectivity tier
       ──> confidence flag from max nearest-neighbour Tanimoto
            (High >= 0.40, Moderate >= 0.30, Low >= 0.20, Outlier < 0.20)
```

- **Reference data:** 794 inhibitors profiled across 464 kinases (141,612 measured activities from ChEMBL + BindingDB).
- **Selectivity metric:** the coverage-robust **S-score** (unlike the Gini coefficient, it is not confounded by how many kinases a compound was assayed against).
- **Applicability domain:** predictions on out-of-domain queries (nearest-neighbour Tanimoto < 0.30) are flagged rather than reported with false confidence.

---

## Repository structure

```
.
├── app/                         # the tool (execution)
│   ├── app.py                   #   Streamlit web interface (single-compound + batch screening)
│   ├── predictor.py             #   core k-NN read-across engine (HuKSARx, morgan_fp)
│   ├── data/                    #   reference data the tool runs on
│   │   ├── compounds_x_kinases_pActivity.csv   # 794 x 464 activity matrix
│   │   ├── kinase_metadata_converted.csv       # kinase gene/family metadata
│   │   └── atlas_background_data.csv            # SMILES, UMAP coords, S-score (for the map view)
│   └── assets/
├── analysis/                    # reproducibility / validation code
│   ├── reproduce_selectivity_analysis.py   # regenerates the selectivity-metric results
│   ├── loco_cv.py               #   leave-one-compound-out cross-validation
│   ├── knn_readacross.py        #   out-of-sample target-recovery (k sweep)
│   └── data/                    #   small inputs for the validation scripts
│       ├── external_panel_30.csv               # 30-compound external test set
│       └── *_huksa_predictions_FROZEN.csv      # frozen, timestamped prediction records
├── screening/
│   └── screen_imppat.py         # large-scale virtual screen of a compound library
├── requirements.txt
├── LICENSE
└── README.md
```

> **Scope note:** this repository contains the **analysis and execution** code only. The figure-,
> table-, and manuscript-generation scripts are intentionally excluded, as are the superseded
> clustering-era experiments.

---

## Installation

```bash
git clone <your-repo-url> huksa && cd huksa
python -m venv .venv && . .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Python 3.11 is recommended (the pinned toolchain was tested on it).

---

## Usage

### 1. Run the web app (execution)
```bash
streamlit run app/app.py
```
Paste or draw a SMILES to get ranked predicted targets with their supporting measured neighbours,
the S-score and selectivity tier, and the confidence flag. A batch mode screens libraries of up to
~100,000 molecules.

### 2. Predict from the command line
```bash
python app/predictor.py "CC1=CC=C(C=C1)..."     # prints a JSON prediction for one SMILES
```

### 3. Reproduce the selectivity-metric analysis
```bash
python analysis/reproduce_selectivity_analysis.py
```
Prints every reported selectivity-metric number (each line maps to a manuscript claim) and writes
CSV outputs to `analysis/repro_outputs/`. Expected headline values include
`r(Gini, coverage) = -0.997`, `S-score coverage r = -0.12`, and `log-entropy median = 0.999`.
The optional Section 6 (K_d validation) is skipped automatically unless you supply the large raw
file (see **Data** below).

### 4. Cross-validation and external validation
```bash
python analysis/loco_cv.py          # leave-one-compound-out ROC-AUC / RMSE over the 794-compound matrix
python analysis/knn_readacross.py   # out-of-sample target recovery (top-1/3/10) and k sweep
```
These were authored to run from the original research tree and reference data by relative path; see
the comments at the top of each script for the exact inputs. The reference matrix they need ships in
`app/data/`, and the frozen prediction records are in `analysis/data/`.

### 5. Virtual screening
```bash
python screening/screen_imppat.py
```
Screens a structure library against the reference set. Point `LIB` (top of the script) at a folder of
`.sdf`/`.mol` files. The IMPPAT natural-product library used in the paper is **not bundled** (it is
~11,700 files); download it from the IMPPAT database if you wish to reproduce that screen.

---

## Data

| File | Bundled? | Notes |
|---|---|---|
| `app/data/compounds_x_kinases_pActivity.csv` | ✅ | 794 × 464 reference activity matrix |
| `app/data/kinase_metadata_converted.csv` | ✅ | kinase gene/family metadata |
| `app/data/atlas_background_data.csv` | ✅ | SMILES + UMAP coords + S-score (map view) |
| `analysis/data/external_panel_30.csv` | ✅ | 30-compound external validation set |
| `analysis/data/*_FROZEN.csv` | ✅ | timestamped frozen prediction records |
| `all_bioactivity_with_censoring.csv` | ❌ | raw long-format bioactivities (~99 MB); only needed for the optional K_d-validation section. Place in `analysis/data/` to enable it. |
| IMPPAT structure library | ❌ | ~11,700 `.sdf`/`.mol` files; download from the IMPPAT database for the virtual screen. |

All bundled bioactivity data derives from the public **ChEMBL** and **BindingDB** databases.

---

## Citation & license

If you use HuKSA, please cite the associated manuscript (in preparation).

**Authors:** Mohammad Umar Saeed and Md Imtaiyaz Hassan.

Released under the **MIT License** — see [LICENSE](LICENSE).
