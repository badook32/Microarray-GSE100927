# Microarray Analysis of GSE100927 for OASL-Related DEG Identification

This repository contains an R script for processing and analyzing the **GSE100927** microarray dataset to identify differentially expressed genes (DEGs) associated with **OASL expression levels** and **arterial site comparisons**. The pipeline includes data loading, normalization, phenotype annotation, DEG identification using `limma`, and probe-to-gene symbol mapping via `biomaRt`.

---
## 🔧 Input Data

* **Platform**: Agilent SurePrint G3 Human GE 8x60K v2
* **Raw Data Format**: `.txt.gz` files for each GSM sample
* **Sample Metadata**: Manually matched using GSM IDs to sample names

---
## 📈 Workflow Overview

### 1. Raw Data Preprocessing

* Load all `.txt.gz` files under the GSE100927 directory
* Extract processed signal values from each file
* Merge individual files into a single expression matrix using `ProbeName` as the key
* Apply log2 transformation to normalize expression values

### 2. Phenotype Annotation

* Assign each sample to an artery site and control status:

  * `Femoral artery`, `Carotid artery`, `Infra-popliteal artery`
  * Corresponding control groups: `*_Control`
* Clean and harmonize group labels for downstream modeling

### 3. Differential Expression Analysis (Limma)

DEGs were computed across the following comparisons using the `limma` pipeline:

* **Femoral artery vs Femoral artery\_Control**
* **Carotid artery vs Carotid artery\_Control**
* **Infra-popliteal artery vs Infra-popliteal artery\_Control**
* **All Disease samples vs All Control samples**

Filtering criteria:

* Adjusted p-value < 0.05 (Benjamini-Hochberg correction)

### 4. Probe Annotation

* Probes were mapped to gene symbols using `biomaRt` with Ensembl v108
* Annotated DEGs were saved for each comparison

### 5. OASL-Based Grouping and DEG Analysis

* Extracted the expression of OASL by probe ID
* Grouped samples into **OASL High** and **OASL Low** based on median expression
* Performed DEG analysis under multiple settings:

  * Across all samples: OASL High vs Low
  * Within **Control** samples only
  * Within **Disease (Other)** samples only

### 6. Visualization

* Boxplot of OASL expression between Control and Disease groups
* Saved as high-resolution TIFF for publication


---
## 📚 References

* GSE100927: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100927)

