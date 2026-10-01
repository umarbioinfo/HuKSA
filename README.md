# HuKSA: Human Kinome Selectivity Atlas

HuKSA is a ligand-based tool for ranking human kinase targets and estimating selectivity from a
query molecule. It uses ECFP4 fingerprints, Tanimoto similarity, and measured kinase activities to
make predictions that can be inspected through their supporting reference compounds.

**Web application:** [huksalab-tool.hf.space](https://huksalab-tool.hf.space)

## How it works

The query is represented as a 2048-bit ECFP4 fingerprint and compared with the reference atlas of
794 compounds and 464 kinases. Target activities are estimated from the five nearest compounds
using a Tanimoto-weighted mean. A target is ranked only when at least two of those neighbours have
a measured activity for it. The measured-panel S-score is the fraction of measured kinases with
pActivity greater than 6.

The maximum query-to-atlas Tanimoto similarity is reported as a structural-support band: High
(at least 0.40), Moderate (at least 0.30), Low (at least 0.20), or Outlier (below 0.20). These are
empirical similarity descriptors, not probabilities of a correct prediction. Missing predictions
are not classified as inactive.

## Validation

The atlas-excluded external benchmark contains 125 kinase inhibitors and probes from KLIFS and
PKIDB. The known target appeared in the top 1 for 24/125 compounds (19.2%), top 3 for 42/125
(33.6%), and top 10 for 69/125 (55.2%). In the High structural-support band, top-10 recovery was
35/42 (83.3%). The support-band results describe this benchmark and should not be interpreted as
per-prediction probabilities.

Internal leave-one-compound-out validation gave pooled ROC-AUC 0.830, precision-recall AUC 0.575,
RMSE 0.647 pActivity units, and Spearman correlation 0.642. Performance declined under scaffold
and similarity exclusions, supporting interpretation as analogue-based read-across rather than
general scaffold hopping. Selectivity prediction and prediction-coverage analyses are available
in the repository outputs.

## Repository contents

```text
.
|-- .github/          GitHub configuration
|-- analysis/         Analysis and validation scripts
|-- app/
|   |-- assets/       Application assets
|   |-- data/         Reference data and map files
|   |-- app.py        Streamlit application
|   `-- predictor.py  Prediction engine and command-line interface
|-- data_manifest/    External-benchmark records and provenance
|-- outputs/          Validation results and tables
|-- screening/        Library-screening code
|-- app.py            Root-level application script
|-- predictor.py      Root-level predictor script
|-- requirements.txt  Pinned Python dependencies
|-- LICENSE
`-- README.md
```

## Run locally

Clone the repository, create and activate a Python virtual environment, then install the pinned
dependencies:

```bash
git clone https://github.com/umarbioinfo/HuKSA.git
cd HuKSA
python -m venv .venv
```

Activate the environment, install the dependencies, and start the app. On Windows PowerShell:

```powershell
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
streamlit run app/app.py
```

On macOS or Linux, use `source .venv/bin/activate` in place of the activation command above.

## Command-line example

Pass a complete SMILES string to the predictor. For example, this runs KU-0063794:

```bash
python app/predictor.py "C[C@@H]1CN(C[C@@H](O1)C)C2=NC3=C(C=CC(=N3)C4=CC(=C(C=C4)OC)CO)C(=N2)N5CCOCC5"
```

The command prints a JSON prediction. The target-ranking portion of the expected output is:

```json
{
  "predicted_targets": [
    {"rank": 1, "gene": "MTOR", "pred_pActivity": 8.7, "neighbours_measured": 4},
    {"rank": 2, "gene": "LRRK2", "pred_pActivity": 6.8, "neighbours_measured": 2},
    {"rank": 3, "gene": "ROS1", "pred_pActivity": 6.7, "neighbours_measured": 2}
  ]
}
```

This is an excerpt of the output; the full JSON also reports structural support and selectivity.
`neighbours_measured` is the number of the five nearest compounds with a measurement for that
kinase.

## Interpretation

- The S-score is conditional on which kinases were measured and on the underlying assay data.
- Structural-support bands summarize similarity to the atlas and are not calibrated probabilities.
- Neighbour activity spread is descriptive and is not a prediction interval.
- External validation is retrospective. HuKSA is intended to support experimental prioritization,
  not replace biochemical kinase profiling.

## Authors and license

Mohammad Umar Saeed, Arunabh Choudhury, Yash Mathur, and Md Imtaiyaz Hassan.

Released under the MIT License. See [LICENSE](LICENSE).
