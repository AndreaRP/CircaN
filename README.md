# CircaN

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R-project](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)

**CircaN** is an R package for robust detection of circadian gene expression using nonlinear modeling.  
It integrates multiple algorithms (CircaN, JTK, MetaCycle) and supports flexible oscillatory patterns, enabling accurate detection in noisy and heterogeneous transcriptomic datasets.

---

## Key Features

- Evaluates **7 distinct oscillatory architectures** (e.g., cosine, triangular) using AIC to ensure optimal fit for heterogeneous data
- Integrates multiple circadian detection methods (**CircaN, JTK, MetaCycle**)  
- Combines statistical evidence using **Fisher’s method** with **Benjamini–Hochberg correction**  
- Handles **noisy, sparse, and heterogeneous datasets**  
- Designed for **reproducible and scalable transcriptomic analysis**  

---

## Installation

```r
# Install from GitHub
install.packages("devtools")
devtools::install_github("AndreaRP/CircaN")
````

---

## Quick Start

```r
library(CircaN)

# Load example data
expression_example <- CircaN::expression_example
rownames(expression_example) <- expression_example$feature
expression_example <- expression_example[,-1]

metadata_example <- CircaN::metadata_example

# Run analysis
results <- full_mode_analysis(
  data = expression_example,
  s2c = metadata_example
)
```

---

## Input Requirements

### Expression Data (`data`)

* Data frame with features as **rownames**
* All columns must be **numeric**
* Each column corresponds to a sample

### Metadata (`s2c`)

Must include:

* `sample`: sample names (matching expression data columns)
* `time`: time point of sample collection
* `ind`: individual/sample identifier

See `metadata_example` for reference.

---

## Parameters 

* `algorithms`: vector of methods to use (`"circan"`, `"jtk"`, `"metacycle"`)
* `circan_mode`: NLS algorithm (`"default"`, `"plinear"`, `"port"`)
* `circan_init_value`: initial period (default: 24)
* `min_per`, `max_per`: period search range (default: 20–28)
* `mc_cycMethod`: algorithms to test in the MetaCycle meta2d (if selected).

---

## Workflow Overview

1. Fit each feature to **7 oscillatory models** using nonlinear regression
2. Select best model using **Akaike Information Criterion (AIC)**
3. Combine statistical evidence across methods using **Fisher’s method**
4. Adjust p-values using **Benjamini–Hochberg correction**

---

## Output

The results object includes:

* Best-fit model per feature
* Estimated **period** and **amplitude**
* Combined **p-values** and adjusted **q-values**

Recommended filtering:

* Adjusted p-value < 0.05
* Goodness-of-fit (R) > 0.7

---

## Use Cases

* Identification of circadian genes in bulk RNA-seq datasets
* Analysis of time-series transcriptomics across tissues or conditions
* Benchmarking circadian detection methods under varying noise conditions
* Integration into larger **multi-omics pipelines**

---

## Performance

CircaN was benchmarked using simulated datasets with varying:

* Noise levels
* Sampling frequencies
* Oscillatory patterns

It demonstrates improved robustness in heterogeneous datasets compared to existing methods.

---

## References & Methodology

CircaN builds upon and integrates established algorithms in the circadian field:

* **JTK_CYCLE**: an efficient non-parametric algorithm for detecting rhythmic components in genome-scale datasets. [Hughes et al. (2010)](https://doi.org/10.1177/0748730410379711) | [Project Website](https://openwetware.org/wiki/HughesLab:JTK_Cycle)
* **MetaCycle**: an integrated R package to evaluate periodicity in large scale data. [Wu et al. (2016)](https://doi.org/10.1093/bioinformatics/btw405) | [Vignette](https://cran.r-project.org/web/packages/MetaCycle/vignettes/implementation.html)

---

## Citation

Rubio-Ponce, A., Ballesteros, I., Quintana, J. A., Solanas, G., Benitah, S. A., Hidalgo, A., & Sánchez-Cabo, F. (2021).
**Combined statistical modeling enables accurate mining of circadian transcription.**
*NAR Genomics and Bioinformatics*, 3(2), lqab031.
[https://doi.org/10.1093/nargab/lqab031](https://doi.org/10.1093/nargab/lqab031)

---

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

---
