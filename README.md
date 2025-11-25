
# Environment & Baseline Notes

## Our Method

This repository contains the official implementation for the paper *“Learning Subgroups with Maximum Treatment Effects without Causal Heuristics”*.  

Our method is built on a simple theoretical insight: **under a partition-based structural causal model, the subgroup with the maximum treatment effect is exactly one of the underlying partitions.**

This allows us to avoid specialized causal heuristics entirely. The algorithm is:

1. Learn a standard CART partition of the feature space (classification or regression).  
2. Estimate subgroup treatment effects honestly on a held-out test split.  
3. Select the leaf with the largest estimated effect.

This pipeline is **theoretically justified, interpretable, and empirically strong**.

```bibtex
@inproceedings{yang2026maximumeffect,
  title={Learning Subgroups with Maximum Treatment Effects without Causal Heuristics},
  author={Yang, Lincen and Li, Zhong and van Leeuwen, Matthijs and Salehkaleybar, Saber},
  booktitle={Proceedings of the AAAI Conference on Artificial Intelligence},
  year={2026}
}
```


## Implementation Summary

* Our method is implemented in **R**.  
* The baselines are a mix of **R** and **Python** code.  
* Because a single environment causes package conflicts, we **highly recommend**:
  * `renv` for managing R packages *per algorithm*.
  * Conda (or another tool) for the Python environment.

---

## How Baselines Are Run

For every baseline we either use the authors’ released source code or an existing R/Python package.  
Each baseline has two main scripts:

```
run_synthetic.R   # or .py
run_semi.R        # or .py
```

---

## R Environment

* R version tested: **4.2.1**  
* OS tested: **macOS** and **Linux**

---

### 1. Semi‑synthetic Dataset Simulator

* `dplyr`
* `data.table`
* `remotes`
* Install simulator:  

  ```r
  remotes::install_github("vdorie/aciccomp/2016")
  ```

### 2. Our Proposed Method

* `rpart`
* `rpart.plot`

### 3. QUINT

* `quint` (+ dependencies)
* `partykit`
* `grid`
* `libcoin`
* `mvtnorm`
* `Formula`  
  Full dependency list: <https://cran.r-project.org/web/packages/quint/>

### 4. Virtual Twins

* `aVirtualTwins` (+ its documented dependencies)  
  <https://cran.r-project.org/web/packages/aVirtualTwins/>

### 5. SIDES

* `rsides` (+ documented dependencies)  
  <https://cran.r-project.org/web/packages/rsides/>

### 6. Causal Tree

* `causalTree` (+ documented dependencies)  
  <https://github.com/susanathey/causalTree>

### 7. Distill Tree

Source code is vendored under `./causalDT-main/` because CRAN install fails on some systems.

Additional packages:

* `rpart`
* `rpart.plot`
* `purrr`
* `R.utils`
* `grf`

Original install instructions:

```r
# install.packages("devtools")
devtools::install_github("tiffanymtang/causalDT", subdir = "causalDT")
```

### 8. Interaction Tree

* Not on CRAN; source released from the authors are publicly at  
  <https://drive.google.com/file/d/1XUIWapeSE3m5VdteTxmyYgSP3A3eeNy4/view>  
* Implementation file: `Functions-IT.R`

---

## Python Baseline — CURLS

We use the authors’ released source.

Add paths at the top of your script if needed:

```python
import sys, os
sys.path.append('./CURLS/osfstorage-archive/')
sys.path.append('./CURLS/osfstorage-archive/pipeline/')
```

* Dependencies list: <https://osf.io/zwp2k/>
* **Note:** CURLS is slow; we ran it on Linux servers with parallel jobs.

Folder structure:

* Main code subfolder: `./curls/osfstorage-archive/pipeline/`
* Original `requirements.txt` in `curls/osfstorage-archive/`
* Our own frozen env: `curls/osfstorage-archive/pipeline/our_requirements.txt`

### Running CURLS

Semi‑synthetic:

```bash
python run_semi.py
# or specify dataset path
python run_semi.py --data_path ../../../datasets/semi_synthetic/ACIC2016_1.csv
```

Synthetic:

```bash
python run_synthetic.py
python run_synthetic.py --data_path ../../../datasets/synthetic/simulate1/n_1000_iter_1.csv
```

