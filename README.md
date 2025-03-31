# RNA-Seq Analysis Project

## Overview
This project analyzes RNA-Seq data to identify tissue-specific and individual-specific gene expression patterns across multiple donors and tissues.

## Project Structure
```
RNA_Seq_Analysis/
├── data/
│   ├── Metadata2025.csv    # Sample metadata
│   └── Counts2025.csv      # Raw count data
├── scripts/
│   └── rna_seq_analysis.R  # Main analysis script
└── results/                # Generated output files
```

## Requirements
- R version 4.0+
- Required R packages:
  - edgeR
  - limma
  - ggplot2
  - DESeq2
  - pheatmap
  - RColorBrewer
  - tidyverse
  - gridExtra
  - stringr

## Analysis Features
1. **Quality Control & Preprocessing**
   - Low expression filtering
   - TMM normalization
   - Sample relationship visualization (PCA & MDS)

2. **Individual Expression Analysis**
   - Individual-specific gene identification
   - Expression pattern classification
   - Cross-individual comparisons

3. **Differential Expression Analysis**
   - Liver vs Brain comparison
   - Multiple FDR thresholds (0.05, 0.01)
   - Volcano plot visualization

## Key Outputs
- `individual_gene_stats.csv`: Individual-specific expression statistics
- `expression_classes_plot.png`: Distribution of expression patterns
- `PCA_plot.png`: Principal Component Analysis visualization
- `MDS_plot.png`: Multi-dimensional Scaling plot
- `liver_vs_brain_DE_genes.csv`: Differential expression results
- `volcano_plot.png`: Differential expression visualization
- `gene_class_boxplot.png`: Expression pattern distributions

## Usage
1. Install required R packages:
```R
install.packages(c("edgeR", "limma", "ggplot2", "DESeq2", "pheatmap", 
                  "RColorBrewer", "tidyverse", "gridExtra", "stringr"))
```

2. Place input files in the `data` directory:
   - `Metadata2025.csv`: Tab-separated metadata file
   - `Counts2025.csv`: Tab-separated count matrix

3. Run the analysis:
```R
source("scripts/rna_seq_analysis.R")
```

4. Results will be generated in the `results` directory

## Expression Classifications
- **Individual-specific**: 4x higher expression in one individual
- **Individual-elevated**: 2x higher expression in one individual
- **Not elevated**: Similar expression across individuals

## Differential Expression Criteria
- FDR < 0.05: Standard threshold
- FDR < 0.01: Stringent threshold
- Log2 fold change used for effect size
