# FAdV-4 Transcriptomic Analysis

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![R](https://img.shields.io/badge/R-4.2+-blue.svg)
![Bioconductor](https://img.shields.io/badge/Bioconductor-3.18+-orange.svg)

> **Integrative transcriptomic analysis reveals conserved host response mechanisms across Fowl Adenovirus Serotype 4 strains**

---

## Overview

This repository contains a complete analysis pipeline for our integrative transcriptomic study of Fowl Adenovirus serotype 4 (FAdV-4) infection. We used **limma combined analysis with batch covariate adjustment** to identify **1,078 differentially expressed genes (DEGs)** from two independent datasets (GSE106839 and GSE299945), with **489 genes validated by cross-dataset consistency assessment**.

### Key Contributions

- Rigorous statistical framework: limma with batch covariate adjustment + cross-dataset validation
- Cross-strain conservation: 489 host response genes conserved across different FAdV-4 strains
- High validation rate: 85.4% direction consistency between datasets
- Therapeutic targets: PDK4, DTL, E2F7, BAMBI, TNFSF10 prioritized
- Fully reproducible: Complete code, data, Docker/Conda environments

---

## Summary of Findings

| Metric | Value |
|--------|-------|
| **Initial common genes** | 10,455 |
| **After filtering** | 8,414 (19.5% removed) |
| **Total DEGs** | 1,078 |
| **Up-regulated** | 459 (42.6%) |
| **Down-regulated** | 619 (57.4%) |
| **Cross-validated DEGs** | 489 |
| **Up-regulated (Cross-validated)** | 165 (33.7%) |
| **Down-regulated (Cross-validated)** | 324 (66.3%) |
| **Direction consistency** | 85.4% (2,732/3,199) |
| **Cross-dataset correlation** | r = 0.536 (logFC) |

---

## Repository Structure

```
FAdV4_LimmaAnalysis/
├── data/
│   ├── processed/
│   │   └── Integrated_Core_DEGs.csv
│   │   └── string_interactions.tsv
│   └── raw/
│       └── source_data/
│           ├── GSE106839/
│           ├── GSE299945/
│           └── sample_metadata.csv
├── results/
│   ├── figures/
│   │   ├── PCA_Comparison.pdf
│   │   ├── Correlation_Heatmap_Raw.pdf
│   │   ├── Correlation_Heatmap_BatchAdjusted.pdf
│   │   ├── Gene_Classification_Summary.pdf
│   │   ├── Enrichment_*.pdf
│   │   ├── STRING_Protein_Interaction.pdf
│   │   ├── WGCNA_*.pdf
│   │   └── ...
│   ├── tables/
│   │   ├── CombinedAnalysis_DEGs.csv (1,078 genes)
│   │   ├── Consistent_DEGs.csv (489 genes)
│   │   ├── CrossDataset_Consistency_Results.csv
│   │   ├── Sample_Information_Full.csv
│   │   ├── WGCNA_Hub_Genes.csv
│   │   ├── WGCNA_Module_Assignment.csv
│   │   ├── WGCNA_Module_Trait_Correlation.csv
│   │   ├── enrichment/
│   │   │   ├── Core_Total_KEGG_Enrichment.csv
│   │   │   ├── Core_Total_BP_GO_Enrichment_BP.csv
│   │   │   └── ...
│   │   └── ppi/
│   ├── intermediate/
│   │   ├── WGCNA_workspace.RData
│   │   ├── Normalized_Expression_BatchAdjusted.rds
│   │   └── ...
│   └── sup_file/
│       ├── Table_S1_AllGenes.csv
│       ├── Table_S2_PPI_Network.csv
│       ├── Table_S3_GO_MF_Enrichment.csv
│       ├── Table_S4_GO_CC_Enrichment.csv
│       └── *.xml (KEGG pathway maps)
├── scripts/
│   ├── config.R
│   ├── Install_Packages.R
│   ├── FAdV-4_Transcriptome_Analysis.R
│   ├── WGCNA_Analysis.R
│   ├── string_interaction_protein.R
│   └── generate_submission_figures.py
├── README.md
├── LICENSE
├── Dockerfile
└── environment.yml
```

---

## Analysis Parameters

### Data Preprocessing
- **Initial common genes**: 10,455
- **Low-expression filter**: Removed 2,041 genes (19.5%)
- **Final genes**: 8,414
- **Normalization**: TMM (edgeR)

### Differential Expression Criteria
- **FDR threshold**: < 0.05 (Benjamini-Hochberg adjusted)
- **logFC threshold**: not applied (continuous values used)

### Batch Effect Correction
- ComBat (sva package) for batch covariate adjustment
- Batch variable: dataset source (GSE106839 vs GSE299945)

---

## Cross-Dataset Consistency

| Dataset | Significant DEGs (FDR<0.05) |
|---------|------------------------------|
| GSE106839 | 4,060 |
| GSE299945 | 6,095 |
| **Overlap** | **3,199** |
| **Direction Consistent** | **2,732 (85.4%)** |

---

## WGCNA Module-Trait Correlations

| Module | Correlation | P-value | Interpretation |
|--------|-------------|---------|----------------|
| ME6 | 0.968 | 2.59×10⁻⁷ | Strong positive |
| ME5 | 0.803 | 1.67×10⁻³ | Moderate positive |
| ME3 | -0.869 | 2.41×10⁻⁴ | Strong negative (contains APOB) |
| ME4 | -0.635 | 2.65×10⁻² | Moderate negative (contains FGA, F13B) |

### WGCNA Parameters
- **Soft threshold power**: 14
- **Scale-free R²**: 0.871
- **Total modules**: 10

---

## Key Pathway Enrichment (Cross-validated DEGs)

### KEGG Pathways

| Pathway | ID | Adjusted P-value |
|---------|-----|------------------|
| Tryptophan metabolism | gga00380 | 8.76×10⁻⁴ |
| Neuroactive ligand-receptor interaction | gga04080 | 1.40×10⁻⁴ |
| Cytokine-cytokine receptor interaction | gga04060 | 1.65×10⁻⁴ |

---

## Top Differentially Expressed Genes

### Up-regulated (Top 5)
| Gene | logFC | Function |
|------|-------|----------|
| CYP1B1 | +5.30 | Cytochrome P450, xenobiotic metabolism |
| PLK3 | +2.32 | Cell cycle kinase |
| CCNO | +1.88 | Cell cycle protein |
| TNFSF10 | +1.75 | Apoptosis signaling |
| DTL | +1.63 | Cell cycle checkpoint |

### Down-regulated (Top 5)
| Gene | logFC | Function |
|------|-------|----------|
| TDO2 | -5.45 | Tryptophan degradation |
| APOB | -5.12 | Lipid transport |
| PPARA | -2.45 | Lipid metabolism |
| PDK4 | -2.68 | Metabolic regulation |
| FGA | -2.31 | Blood coagulation |

---

## PPI Network Analysis

- **Network nodes**: 489 genes (Cross-validated DEGs)
- **Interaction score threshold**: > 0.7
- **Hub genes**: Top 10% by degree (n=30)
- **APOB degree**: 94 (highest connectivity)
- **FGA degree**: 84
- **F13B degree**: 76

---

## Druggability Analysis

- 30 PPI hub genes analyzed
- 21 genes (70.0%) with DisGeNET disease annotations
- 18/21 (85.7%) belong to druggable genome
- **APOB**: Only clinically actionable gene identified

---

## Dependencies

### R Packages (>= 4.2)
- limma, edgeR, DESeq2
- sva, batchelor
- clusterProfiler, DOSE, enrichplot
- WGCNA, igraph
- ggplot2, pheatmap

### Python
- pandas, matplotlib, seaborn

### Docker/Conda
- Full environment specification in Dockerfile and environment.yml

---

## Reproducibility

All analyses can be reproduced using:
1. **Docker**: `docker build -t fadv4-analysis .`
2. **Conda**: `conda env create -f environment.yml`

---

## License

MIT License

---

## Citation

[To be added upon publication]