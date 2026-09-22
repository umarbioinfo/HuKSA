# HuKSA - Human Kinome Selectivity Atlas

HuKSA is a transparent ligand-based framework for prioritizing human kinase targets and estimating
measured-kinome selectivity from a SMILES string. It uses molecular fingerprint similarity and
Tanimoto-weighted k-nearest-neighbour read-across over measured kinase activities. Each prediction
is accompanied by the reference compounds and measurements that support it.

**Web application:** https://huksalab-tool.hf.space

## Method

```text
SMILES -> 2048-bit ECFP4 fingerprint -> Tanimoto similarity to 794 atlas compounds
       -> similarity-weighted k-nearest-neighbour read-across with k = 5
       -> predicted pActivity and ranked kinase targets
       -> measured-panel S-score and selectivity tier
       -> empirical structural-support band from maximum atlas similarity
```

- **Reference atlas:** 794 compounds, 464 kinases, and 141,612 measured activities.
- **Target prediction:** a Tanimoto-weighted mean is calculated from the five nearest atlas
  compounds. A kinase is rankable only when at least two neighbours have a measured activity for
  that kinase.
- **Selectivity:** the S-score is the fraction of measured kinases with pActivity greater than 6.
  It avoids the direct assay-count artefact caused by padding unmeasured entries, but remains
  conditional on panel composition and assay type.
- **Structural support:** High is at least 0.40, Moderate at least 0.30, Low at least 0.20, and
  Outlier below 0.20 maximum atlas Tanimoto similarity. These bands describe analogue support and
  are not calibrated probabilities.
- **Evidence display:** the application reports measured-neighbour identities, support counts,
  weighted neighbour activity spread, and the number of rankable kinases for each query.

## Validation

### External target recovery

The main external benchmark contains 125 atlas-excluded KLIFS and PKIDB kinase inhibitors and
probes. Structural overlap checks included standardized structures, full InChIKeys, InChIKey
skeleton blocks, canonical tautomers, Murcko scaffolds, and fingerprint similarity.

| Metric | Result |
|---|---:|
| Top-1 target recovery | 24/125, 19.2% |
| Top-3 target recovery | 42/125, 33.6% |
| Top-10 target recovery | 69/125, 55.2% |
| High-support top-10 recovery | 35/42, 83.3% |

The structural-support bands are empirical descriptors. They do not give the probability that an
individual target prediction is correct.

### Internal target validation

Standard leave-one-compound-out validation produced pooled ROC-AUC 0.830, precision-recall AUC
0.575, RMSE 0.647 pActivity units, and Spearman correlation 0.642. ROC-AUC decreased to 0.807 after
Murcko-scaffold exclusion and to 0.775, 0.744, and 0.713 after excluding neighbours at Tanimoto
similarity thresholds of 0.50, 0.40, and 0.30. HuKSA should therefore be interpreted as an
analogue-supported read-across method rather than a general scaffold-hopping model.

### Selectivity validation

Direct prediction of the pActivity 6 S-score gave RMSE 0.123, Spearman correlation 0.562, and tier
accuracy 55.0%. Under Murcko-scaffold exclusion, the corresponding values were 0.134, 0.499, and
51.8%.

Among 74 compounds measured on a common 360-kinase panel, full-panel and common-panel S-scores gave
Spearman correlation 0.998 and tier agreement 98.6%. Repeated matched 100-kinase panels gave median
tier agreement of 87.8% for random panels and 86.5% after kinase-family balancing.

### Prediction coverage

Of 368,416 compound-kinase combinations in the internal grid, 158,656 were rankable under the
two-neighbour rule, giving 43.1% coverage. Median per-compound coverage was 30.6% and median
per-kinase coverage was 35.0%. The median weighted neighbour standard deviation was 0.239 pActivity
units. This value describes local disagreement among measured neighbours and is not a calibrated
prediction interval.

## Repository contents

```text
.
|-- app/
|   |-- app.py                 Streamlit interface
|   |-- predictor.py           k = 5 read-across and evidence engine
|   |-- assets/
|   `-- data/                  activity matrix, kinase metadata, and map data
|-- analysis/                  validation and reproducibility scripts
|-- data_manifest/             external-set protocol, records, and checksum
|-- outputs/                   validation tables and summaries
|-- figures/                   manuscript and supplementary figures
|-- supplementary_data/        per-compound and per-kinase coverage tables
|-- tests/                     regression and output-integrity tests
|-- requirements.txt
|-- LICENSE
`-- README.md
```

## Installation

```bash
git clone https://github.com/umarbioinfo/HuKSA.git
cd HuKSA
python -m venv .venv
```

Activate the environment and install the dependencies.

```bash
python -m pip install -r requirements.txt
```

## Usage

Run the web application locally.

```bash
streamlit run app/app.py
```

Run a command-line prediction.

```bash
python app/predictor.py "CC1=CC=C(C=C1)..."
```

Run the test suite.

```bash
python -m pytest tests
```

## Interpretation limits

- IC50, Ki, and Kd measurements are not assumed to be interchangeable.
- Missing kinase predictions are omitted from rankings and are not classified as inactive.
- Structural-support bands are empirical similarity descriptors, not probabilities.
- Neighbour standard deviation is a local activity-spread measure, not a prediction interval.
- External validation is retrospective rather than temporal or prospective.
- The natural-product screen is a computational prioritization exercise, not biological validation.
- HuKSA supports experimental prioritization and does not replace biochemical kinase profiling.

## Authors and license

Mohammad Umar Saeed, Arunabh Choudhury, Yash Mathur, and Md Imtaiyaz Hassan.

Released under the MIT License. See [LICENSE](LICENSE).
