# Deep Learning for Biology — TCGA-BRCA Transcriptomic Analysis & Predictive Modeling

Supervised Multi-Layer Perceptron (MLP) project developed for the **UM5BIP12 — AI & Deep Learning for Biology** teaching unit (Master 2 BMC — Artificial Intelligence & Data Analysis, Sorbonne Université).

---

## 1. Dataset Overview & Provenance

All data tables are sourced from the **GDC TCGA Breast Cancer (`TCGA-BRCA`)** cohort hosted on the [UCSC Xena Data Hub](https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Breast%20Cancer%20(BRCA)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443).

| Modality / Asset | Xena Dataset ID | GDC Release & Version | Dimensions / Key Details | Link |
| :--- | :--- | :--- | :--- | :--- |
| **Gene Expression** | `TCGA-BRCA.star_counts.tsv` | GDC Rel. 40.0 (`2024-05-20`) | **60,660 loci × 1,226 samples**; Illumina RNA-seq STAR counts, aliquot-averaged, $\log_2(\text{count} + 1)$ transformed | [Xena Link](https://xenabrowser.net/datapages/?dataset=TCGA-BRCA.star_counts.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) |
| **Gene Annotation** | `gencode.v36.annotation.gtf.gene.probemap` | GENCODE v36 (`2023-08-07`) | Human genome assembly **`hg38` (`GRCh38.p13`)**; maps Ensembl Gene IDs (`ENSG...`) to HGNC gene symbols | [GTF Source](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz) |
| **Clinical Phenotype** | `TCGA-BRCA.clinical.tsv` | GDC Rel. 41.0 (`2024-09-07`) | Patient demographics, tumor staging, receptor status, and molecular annotations | [Xena Link](https://xenabrowser.net/datapages/?dataset=TCGA-BRCA.clinical.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) |
| **Survival Outcomes** | `TCGA-BRCA.survival.tsv` | GDC Rel. 40.0 (`2024-08-03`) | Overall survival (`OS`) and progression-free interval (`PFI`) event indicators and follow-up times | [Xena Link](https://xenabrowser.net/datapages/?dataset=TCGA-BRCA.survival.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) |

---

## 2. Exploratory Data Analysis & Preprocessing (`notebooks/EDA.ipynb`)

### 2.1 Cohort Composition & Library-Size Normalization
* **Sample Breakdown:** The expression matrix spans **1,226 biopsies** stratified by TCGA barcode tissue codes into **Primary Tumors (`-01`)**, **Solid Tissue Normals (`-11`)**, and **Metastatic samples (`-06`, $n = 7$)**.
* **Sequencing Depth Correction:** Reversing Xena's $\log_2(\text{count} + 1)$ transformation revealed a **5.93× variation** in raw library size across biopsies (**19.2M** to **114.0M** mapped reads). Applying sample-level Counts Per Million normalization ($\log_2(\text{CPM} + 1)$) aligned sample expression medians and eliminated sequencing-depth confounding from downstream axes ($r = -0.040$ with PC1, $r = 0.019$ with PC2).

### 2.2 Three-Tier Mean–Variance Transcriptome Filtering
To mitigate the $p \gg n$ dimensionality regime prior to exploratory analysis, all **60,660** GENCODE v36 loci were partitioned across a two-stage filtering funnel:
1. **Low-Expression Noise (`43,464` loci — 71.7%):** Loci failing a detection threshold of $\log_2(\text{CPM} + 1) > 1$ in at least **20%** of biopsies, removing unexpressed background noise and rare-subgroup outliers.
2. **Stable Expressed / Housekeeping Genes (`15,196` genes — 25.1%):** Ubiquitously expressed genes falling below the variance cutoff ($\text{var} < 1.54$), excluded prior to unit-variance scaling (`StandardScaler`) so constitutive genes do not dilute tumor heterogeneity.
3. **Top Highly Variable Genes (`2,000` HVGs — 3.3%):** Highest-variance expressed genes ($\text{var} \ge 1.54$) retained for exploratory dimensionality reduction.

### 2.3 Exploratory PCA & Biological Validation
Standardized PCA across the top 2,000 HVGs concentrates **27.6%** of total variance in the first two components (**PC1: 14.7%**, **PC2: 12.9%**) and **58.5%** across the top 20 components, capturing two primary biological axes without institutional batch distortion:
* **PC1 (14.7%) — Normal vs. Malignant Architecture:** Cleanly separates Solid Tissue Normal (`-11`) biopsies—driven by myoepithelial/stromal markers (*DMD*, *SFRP1*, *ANXA1*, *CAV2*, *CRYAB*)—from Primary (`-01`) and Metastatic (`-06`) tumors (*COX6C*, *SLC9A3R1*, *STARD10*, *BCAS4*).
* **PC2 (12.9%) — Proliferative/Basal vs. Luminal Subtype Split:** Separates a rapidly proliferating tumor sub-cluster driven by mitotic cell-cycle regulators (*CCNE1*, *CDC20*, *AURKB*, *CDCA8*, *KIF2C*, *NDC80*, *MCM10*, *CDCA3*) from hormone-receptor-regulated Luminal tumors (*ERBB4*, *AR*, *SCUBE2*, *IL6ST*, *CCDC170*, *TBC1D9*).
* **Hospital Tissue Source Site (`TSS`) Audit:** The six largest collection centers (**BH**: 227, **E2**: 104, **A2**: 102, **A8**: 84, **D8**: 78, **E9**: 77; 54.8% of total samples) intermingle homogeneously across tumor subtypes, confirming the absence of strong hospital-specific batch clusters (with `BH` enrichment in the normal cluster reflecting paired tumor/normal sampling).

---

## 3. Repository Structure

```text
├── data/                   # Local UCSC Xena raw files (git-ignored)
├── notebooks/
│   └── EDA.ipynb           # Transcriptomic & clinical exploratory data analysis
├── src/                    # Preprocessing, dataset loaders, baselines, and MLP modules
└── README.md
```

---

## 4. Roadmap

- [x] **Phase 1 — Transcriptomic EDA:** Library-size normalization, HVG filtering funnel, PCA biological validation, and TSS batch audit.
- [ ] **Phase 2 — Clinical & Survival Audit:** Align `TCGA-BRCA.clinical.tsv` and `TCGA-BRCA.survival.tsv`, audit missingness/class balance, and finalize the prediction target ($y$).
- [ ] **Phase 3 — Leakage-Free Splitting:** Implement patient-level grouped train/validation/test splits (`TCGA-XX-XXXX` participant IDs) with inside-fold feature selection and scaling.
- [ ] **Phase 4 — Modeling & Evaluation:** Benchmark regularized linear/tree baselines against a regularized PyTorch MLP (`Dropout`, `Weight Decay`, `Early Stopping`) and interpret learned gene attributions.