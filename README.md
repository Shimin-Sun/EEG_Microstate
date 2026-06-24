
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
**Script:** `ami_analysis.R`

Computes pairwise Adjusted Mutual Information (AMI) between microstate label sequences across six frequency bands: *Broadband, Delta, Theta, Alpha, Beta, Gamma*.  
Statistical significance is assessed via **permutation testing** (1,000 shuffles) with **FDR correction** for multiple comparisons.

### 2. Co‑occurrence Rate (ΔP)
**Script:** `cooccurrence_analysis.R`

Calculates row‑wise and column‑wise conditional probability differences (ΔP) for each broadband–narrowband microstate pair.  
FDR correction is applied independently for each frequency‑band combination.

### 3. Transfer Entropy (TE)
**Script:** `transfer_entropy.R`

Estimates effective transfer entropy between broadband and narrowband microstate label sequences, quantifying directed information flow.  
Statistical significance is determined via **bootstrap resampling** (1,000 iterations).

### 4. Machine Learning Classification
**Script:** `classification_pipeline.py`

Evaluates the discriminative power of cross‑frequency coupling features using a **nested cross‑validation** scheme (5 outer × 3 inner folds).

- **Feature selection methods:** RFE, MRMR, ReliefF, Lasso, Random Forest  
- **Classifier:** Random Forest  
- **Performance metrics:** Balanced accuracy, sensitivity, specificity, AUC

The pipeline automatically handles missing data (median imputation) and scaling (`RobustScaler`). Results include feature selection frequency and mean importance scores.

## Usage

### Prerequisites
- **R** (≥ 4.0) with packages: `aricode`, `infotheo`, `RTransferEntropy`, etc.  
- **Python** (≥ 3.8) with packages: `pandas`, `numpy`, `scikit‑learn`, `skrebate`, `openpyxl`

Install Python dependencies with:
```bash
pip install pandas numpy scikit-learn skrebate openpyxl
