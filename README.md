
---

## Requirements

### R
- dplyr
- reshape2
- infotheo
- aricode
- RTransferEntropy
- openxlsx
- parallel
- doParallel
- foreach

### Python
- Python ≥ 3.8
- numpy
- pandas
- scikit‑learn
- matplotlib
- seaborn

---

## Analysis Steps

### 1. Adjusted Mutual Information (AMI)
`ami_analysis.R`

Computes pairwise AMI between microstate label sequences across six frequency bands (Broadband, Delta, Theta, Alpha, Beta, Gamma). Statistical significance is assessed via permutation testing (1,000 shuffles) with FDR correction.

### 2. Co‑occurrence Rate (ΔP)
`cooccurrence_analysis.R`

Calculates row‑wise and column‑wise conditional probability differences (ΔP) for broadband‑narrowband microstate pairs. FDR correction applied per frequency‑band combination.

### 3. Transfer Entropy (TE)
`transfer_entropy.R`

Computes effective transfer entropy between broadband and narrowband microstate label sequences to quantify directed information flow. Statistical significance assessed via bootstrap resampling (1,000 iterations).

### 4. Machine Learning Classification
`classification_pipeline.py`

Evaluates discriminative performance using nested cross‑validation (5 outer × 3 inner folds):
- Feature selection: RFE, MRMR, ReliefF, Lasso, Random Forest
- Classifiers: Random Forest
- Metrics: Balanced accuracy, sensitivity, specificity, AUC

---

## Usage

Update the file paths in each script to your local data directory:

```r
base_path <- "path/to/EC_sequences"
output_dir <- "path/to/output"
