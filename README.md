# FAdV-4 Transcriptomic Analysis

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![R](https://img.shields.io/badge/R-4.2+-blue.svg)
![Bioconductor](https://img.shields.io/badge/Bioconductor-3.18+-orange.svg)

> **Integrative transcriptomic analysis reveals conserved host response mechanisms across Fowl Adenovirus Serotype 4 strains**

---

## 📖 Overview

This repository contains a complete analysis pipeline for our integrative transcriptomic study of Fowl Adenovirus serotype 4 (FAdV-4) infection. We used **limma combined analysis with batch covariate adjustment** to identify **540 differentially expressed genes (DEGs)** from two independent datasets (GSE106839 and GSE299945), with **489 genes (90.6%)** validated by cross-dataset consistency assessment.

### 🎯 Key Contributions

- ✅ **Rigorous statistical framework**: limma with batch covariate adjustment + cross-dataset validation
- ✅ **Cross-strain conservation**: 489 host response genes conserved across different FAdV-4 strains
- ✅ **High validation rate**: 90.6% overlap between combined analysis and cross-dataset consistent genes
- ✅ **Therapeutic targets**: PDK4, DTL, E2F7, BAMBI, TNFSF10 prioritized
- ✅ **Fully reproducible**: Complete code, data, Docker/Conda environments

---

## 📊 Summary of Findings

| Metric | Value |
|--------|-------|
| **Initial common genes** | 10,455 |
| **After filtering** | 8,414 (19.5% removed) |
| **Part2 DEGs (Combined)** | 540 |
| **Part3 Consistent DEGs** | 489 (90.6% of Part2) |
| **Up-regulated (Consistent)** | 179 genes |
| **Down-regulated (Consistent)** | 310 genes |
| **Direction consistency** | 100% (489/489) |
| **Cross-dataset correlation** | r = 0.550 (logFC) |

---

## 🗂️ Repository Structure

```
FAdV4_LimmaAnalysis/
├── data/
│   ├── processed/              # Processed data files
│   │   └── string_interactions.tsv
│   └── raw/                    # Raw data (not in Git)
│       └── source_data/
│           ├── GSE106839/
│           ├── GSE299945/
│           └── sample_metadata.csv
├── results/
│   ├── figures/               # Publication-quality PDFs
│   │   ├── 1_PCA_Comparison.pdf
│   │   ├── 2_Correlation_Heatmap_Raw.pdf
│   │   ├── 3_Correlation_Heatmap_BatchAdjusted.pdf
│   │   ├── Gene_Classification_Summary.pdf
│   │   ├── Enrichment_*.pdf
│   │   ├── STRING_Protein_Interaction.pdf
│   │   └── ...
│   ├── tables/                # Analysis result CSVs
│   │   ├── CombinedAnalysis_DEGs.csv (540 genes)
│   │   ├── Consistent_DEGs.csv (489 genes)
│   │   ├── CrossDataset_Consistency_Results.csv
│   │   ├── STRING_Protein_Interaction.csv
│   │   ├── Sample_Information_Full.csv
│   │   ├── enrichment/        # GO/KEGG enrichment results
│   │   ├── ppi/              # Protein-protein interaction data
│   │   └── ...
│   ├── intermediate/          # Intermediate RDS files
│   └── sup_file/             # Supplementary KEGG pathway maps
├── scripts/
│   ├── config.R              # Configuration file
│   ├── Install_Packages.R    # Package installation script
│   ├── Part1_DataPreprocessing_QC.R
│   ├── Part2_Combined_Analysis_Integrated.R
│   ├── Part3_CrossDataset_Consistency.R
│   ├── Part4_Enrichment_Analysis_Unified.R
│   ├── Part5_Complete_Visualization.R
│   ├── Part6_Sensitivity_Analysis.R
│   └── string_interaction_protein.R
├── README.md
├── LICENSE
├── CODE_REVIEW_REPORT.md     # Detailed code review report
├── Dockerfile
└── environment.yml
```

---

## 🔬 Analysis Parameters

### Data Preprocessing
- **Initial common genes**: 10,455
- **Low-expression filter**: Removed 2,041 genes (19.5%)
- **Final genes**: 8,414
- **Normalization**: TMM (edgeR)

### Differential Expression Criteria
- **FDR threshold**: < 0.05 (Benjamini-Hochberg adjusted)
- **Log2 fold-change**: |logFC| > 1.5 (based on power analysis, ~70-75% power)
- **Cross-dataset consistency**: Required (same direction in both datasets)
- **Validation rate**: 90.6% (489/540 genes validated)

### Statistical Methods
- **Primary analysis (Part2)**:
  - limma combined analysis with batch covariate adjustment
  - Model: ~0 + dataset_factor + group_factor
  - Contrast: group_factorFAdV4 vs Control (adjusted for batch)
- **Cross-dataset validation (Part3)**:
  - Independent analysis of each dataset
  - Intersection of significant genes (FDR < 0.05 in both)
  - Direction consistency check
- **Multiple testing correction**: Benjamini-Hochberg (FDR)

### Statistical Power
- **Sample size**: n = 3 per group per dataset (combined n = 6 per group)
- **Power analysis (Part1)**:
  - n = 3 provides ~70-75% power at |logFC| > 1.5
  - 80% power requires |logFC| ≥ 1.75
- **Threshold justification**:
  - logFC > 1.5 balances statistical rigor with sensitivity
  - Higher than logFC > 1.0 (~60% power)
  - Supported by sensitivity analysis (Part6)

---

## 🧬 Top Consistent Differentially Expressed Genes

| Gene | avg_logFC | Regulation | Function |
|------|-----------|------------|----------|
| **APOB** | -4.09 | Down | Lipid metabolism |
| **KNG1** | -3.91 | Down | Inflammation regulation |
| **CYP1B1** | +4.29 | Up | Xenobiotic metabolism |
| **TDO2** | -4.04 | Down | Tryptophan metabolism |
| **EOMES** | +4.61 | Up | Transcription factor |
| **DTL** | -1.51 | Down | DNA replication licensing |
| **E2F7** | -1.71 | Down | Transcription factor |
| **PDK4** | -2.68 | Down | Glycolysis regulation |
| **BAMBI** | +3.32 | Up | TGF-β signaling |
| **RGS7BP** | +3.00 | Up | Cell signaling |

> Note: avg_logFC represents the average fold-change across both datasets. All genes listed have FDR < 0.05 in both datasets and 100% direction consistency.

---

## 🚀 Installation

### Option 1: Using Docker (Recommended)

```bash
# Build Docker image
docker build -t fadv4-analysis:latest .

# Run analysis (mount data and results directories)
docker run -v $(pwd)/data:/app/data \
           -v $(pwd)/results:/app/results \
           fadv4-analysis:latest \
           Rscript /app/scripts/Part1_DataPreprocessing_QC.R
```

### Option 2: Using Conda

```bash
# Create conda environment
conda env create -f environment.yml
conda activate fadv4-meta

# Verify installation
R --version
```

### Option 3: Manual Installation

```r
# Install R packages
install.packages(c("tidyverse", "data.table", "ggplot2", "pheatmap"))

# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "clusterProfiler",
                       "org.Gg.eg.db", "DOSE", "enrichplot"))
```

---

## 💻 Usage

### Quick Start

```bash
# Run complete pipeline sequentially
cd /path/to/FAdV4_LimmaAnalysis

Rscript scripts/Part1_DataPreprocessing_QC.R
Rscript scripts/Part2_Combined_Analysis_Integrated.R
Rscript scripts/Part3_CrossDataset_Consistency.R
Rscript scripts/Part4_Enrichment_Analysis_Unified.R
Rscript scripts/Part5_Complete_Visualization.R
Rscript scripts/Part6_Sensitivity_Analysis.R
Rscript scripts/string_interaction_protein.R
```

### Analysis Pipeline

#### Step 1: Data Preprocessing & QC
```r
source("scripts/Part1_DataPreprocessing_QC.R")

# Outputs:
# - PCA plots (before/after batch correction)
# - Correlation heatmaps
# - Density distributions
# - Filtered expression matrix (8,414 genes)
# - Power analysis report
```

#### Step 2: Combined Analysis with Batch Covariate
```r
source("scripts/Part2_Combined_Analysis_Integrated.R")

# Parameters:
FDR_THRESHOLD <- 0.05
LOGFC_THRESHOLD <- 1.5

# Outputs:
# - CombinedAnalysis_DEGs.csv (540 genes)
# - All genes with statistics
```

#### Step 3: Cross-Dataset Consistency Assessment
```r
source("scripts/Part3_CrossDataset_Consistency.R")

# Parameters (same as Part2 for consistency):
FDR_THRESHOLD <- 0.05
LOGFC_THRESHOLD <- 1.5

# Outputs:
# - Consistent_DEGs.csv (489 genes)
# - CrossDataset_Consistency_Results.csv
# - Consistency metrics and visualizations
```

#### Step 4: Enrichment Analysis
```r
source("scripts/Part4_Enrichment_Analysis_Unified.R")

# Outputs:
# - GO enrichment results (BP, MF, CC)
# - KEGG pathway analysis
# - STRING network input files
```

#### Step 5: Complete Visualization
```r
source("scripts/Part5_Complete_Visualization.R")

# Outputs:
# - Volcano plots
# - Gene classification summary
# - Heatmaps
# - Box plots
# - Venn diagrams
# - KEGG pathway maps
```

#### Step 6: Sensitivity Analysis
```r
source("scripts/Part6_Sensitivity_Analysis.R")

# Outputs:
# - Sensitivity analysis plots
# - Robustness assessment across multiple FDR/logFC thresholds
```

#### Step 7: STRING Protein Interaction Network
```r
source("scripts/string_interaction_protein.R")

# Outputs:
# - STRING_Protein_Interaction.pdf
# - STRING_Protein_Interaction.csv
```

---

## 📈 Results Summary

### Core Findings

#### 1. Cross-Strain Conserved Host Response
- **540 DEGs** identified using limma combined analysis (Part2)
- **489 genes (90.6%)** validated by cross-dataset consistency assessment (Part3)
- **100% direction consistency** for validated genes
- **Correlation**: r = 0.550 between dataset logFC values

#### 2. Metabolic Reprogramming
- **Key pathways**: Lipid metabolism, glycolysis
- **Key genes**: APOB (-4.09), KNG1 (-3.91), PDK4 (-2.68)
- **Mechanism**: Metabolic shift to support viral replication

#### 3. Cell Cycle Dysregulation
- **Key genes**: DTL, E2F7, CCNE2
- **Pathway**: G1/S transition checkpoint
- **Mechanism**: FAdV-4 hijacks host cell cycle machinery

#### 4. Immune Modulation
- **Key genes**: BAMBI (+3.32), TNFSF10
- **Pathway**: TGF-β signaling, apoptosis regulation
- **Mechanism**: Host immune response modulation

---

## 📦 Output Files

### Figures

**Quality Control (Part1):**
- 1_PCA_Comparison.pdf
- 2_Correlation_Heatmap_Raw.pdf
- 3_Correlation_Heatmap_BatchAdjusted.pdf
- 4_Density_Distributions.pdf
- 5_Boxplots.pdf
- Statistical_Power_Analysis.pdf

**Analysis Results (Part5):**
- Gene_Classification_Summary.pdf
- Enrichment_GO_BP_Top10.pdf
- Enrichment_KEGG_Top10.pdf
- Figure_Heatmap_Top50_Clustered.pdf
- Figure_Boxplot_Top20_DEGs.pdf
- Figure_VennDiagram.pdf

**Network Analysis:**
- STRING_Protein_Interaction.pdf

**Sensitivity Analysis (Part6):**
- Sensitivity_Analysis_Plots.pdf

**KEGG Pathway Maps (results/sup_file/):**
- gga00380.FAdV4.png (Tryptophan metabolism)
- gga04080.FAdV4.png (Oocyte meiosis)
- gga04060.FAdV4.png (Cytokine-cytokine receptor interaction)
- gga03320.FAdV4.png (PPAR signaling pathway)
- gga04920.FAdV4.png (Adipocytokine signaling pathway)

### Tables

**DEG Lists:**
- CombinedAnalysis_DEGs.csv (540 genes from Part2)
- Consistent_DEGs.csv (489 genes from Part3)
- CombinedAnalysis_AllGenes.csv (8,414 genes with full statistics)
- CrossDataset_Consistency_Results.csv (8,414 genes with cross-dataset statistics)

**Enrichment Results (results/tables/enrichment/):**
- Core_Total_GO_Enrichment_BP.csv
- Core_Total_GO_Enrichment_MF.csv
- Core_Total_GO_Enrichment_CC.csv
- Core_Total_KEGG_Enrichment.csv

**Network Data:**
- STRING_Protein_Interaction.csv (280 nodes)
- ppi/STRING_Input_PrimaryDEGs.csv

**Sample Information:**
- Sample_Information_Full.csv (with correct grouping)
- Sample_Information_Raw.csv

**Statistical Analysis:**
- Table1_Top30_DEGs.csv
- Table2_KEGG_Pathways.csv
- Table3_CrossDataset_Summary.csv

**Sensitivity Analysis:**
- Sensitivity_Analysis_Results.csv

**Intermediate Files:**
- Normalized_Expression_For_Statistics.rds (for DEG analysis)
- Normalized_Expression_BatchAdjusted.rds (for visualization)

---

## 🎯 Therapeutic Targets

### Prioritized for Further Study

| Target | avg_logFC | Rationale |
|--------|-----------|-----------|
| **PDK4** | -2.68 | Metabolic reprogramming; existing inhibitors |
| **DTL** | -1.51 | DNA replication; essential for viral proliferation |
| **E2F7** | -1.71 | Transcription factor; master cell cycle regulator |
| **BAMBI** | +3.32 | Immune modulation; TGF-β signaling |
| **TNFSF10** | -1.53 | Apoptosis modulation; soluble factor |
| **APOB** | -4.09 | Lipid metabolism; highly down-regulated |
| **CYP1B1** | +4.29 | Xenobiotic metabolism; highly up-regulated |

---

## 📚 Citation

If you use this code, data, or findings in your research, please cite:

```bibtex
@article{fadv4_transcriptomics_2025,
  title={Integrative transcriptomic analysis reveals conserved host response mechanisms across Fowl Adenovirus Serotype 4 strains},
  author={[Your Name] and [Co-authors]},
  journal={Briefings in Bioinformatics},
  year={2025},
  publisher={Oxford University Press},
  note={Integrative analysis using limma with batch covariate adjustment}
}
```

---

## 🤝 Contributing

We welcome contributions! Areas for improvement:

- **Additional FAdV-4 datasets** - Expand meta-analysis
- **Time-course analysis** - Dynamic host responses
- **Proteomics integration** - Multi-omics validation
- **Cross-species comparison** - Compare with other adenoviruses

Please ensure:
1. Code follows existing style conventions
2. All scripts have proper documentation
3. Results are reproducible
4. Pull requests are made to the `main` branch

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 📧 Contact

- **GitHub**: https://github.com/Tygg199456/Across-FAdV-4-strains-Analysis

---

## 🙏 Acknowledgments

### Data Sources
- **GEO Datasets**: GSE106839, GSE299945
- **PPI Database**: STRING v11.5 (Gallus gallus)
- **Genome Annotation**: Ensembl (Gallus gallus)

### Software & Tools
- **R/Bioconductor**: limma, edgeR, clusterProfiler, org.Gg.eg.db
- **Network Analysis**: STRING database
- **Visualization**: ggplot2, pheatmap, ggrepel

### Statistical Methods
- **limma**: Ritchie et al. (2015) Nucleic Acids Research 43(7):e47
- **TMM normalization**: Robinson & Oshlack (2010) Genome Biology 11:R25
- **BH correction**: Benjamini & Hochberg (1995)

---

## ⚠️ Limitations

1. **Sample size**: n=3 per group per dataset limits statistical power for moderate effect sizes
2. **Dataset availability**: Only two public datasets available
3. **Species-specific data**: Gallus gallus annotations less complete than model organisms
4. **Cross-sectional design**: Single time-point (12hpi) limits dynamic insights
5. **Computational predictions**: Network analysis requires experimental confirmation

---

## 📈 Changelog

### v2.0.0 (2026-02-10) - Current Version
- **Fixed sample grouping**: Corrected metadata grouping for both datasets
  - GSE106839: samples 1-3 are FAdV4, samples 4-6 are Control
  - GSE299945: samples 1-3 are Control, samples 4-6 are FAdV4
- **Updated results**:
  - Part2 (Combined): 540 DEGs (FDR < 0.05, |logFC| > 1.5)
  - Part3 (Consistent): 489 DEGs (90.6% validation rate)
- **Fixed hardcoded template**: Part1_DataPreprocessing_QC.R no longer has incorrect hardcoded grouping
- **Code review**: Comprehensive code review and validation performed
- **Updated documentation**: README, Dockerfile, environment.yml updated

### v1.1.0 (2025-02-07)
- **Updated logFC threshold**: 1.0 → 1.5 (based on power analysis)
- **Clarified methods**: limma combined analysis (not meta-analysis)
- **Enhanced documentation**: Updated all parameter references

### v1.0.0 (2025-02-05)
- Initial release
- Complete analysis pipeline
- Power analysis and sensitivity analysis
- Cross-dataset consistency assessment
- GO/KEGG enrichment analysis
- Docker and Conda environments

---

## 🔗 Links

- **GitHub**: https://github.com/Tygg199456/Across-FAdV-4-strains-Analysis
- **Issues**: https://github.com/Tygg199456/Across-FAdV-4-strains-Analysis/issues
- **GEO Datasets**: https://www.ncbi.nlm.nih.gov/geo/

---

<div align="center">

**If you find this work useful, please consider giving it a ⭐️!**

</div>
