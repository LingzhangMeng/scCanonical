# scCanonical: An R Package for Canonical Marker Identification in scRNA-seq Analysis

## Overview
The `scCanonical` R package is designed to identify and rank canonical markers for downstream cluster annotations in single-cell RNA sequencing (scRNA-seq) analysis. It is specifically developed to work with integrated Seurat objects created using the `Seurat` R package. The package facilitates the calculation of conserved markers across conditions, computes specificity scores, and visualizes results to aid in the interpretation of scRNA-seq data.

## Why you join before computing specificity
FindConservedMarkers() produces a meta-analysis table (p-values across groups), but often doesn’t include avg_log2FC, pct.1, pct.2. FindAllMarkers() does include those columns (computed broadly cluster vs. rest). Joining on (gene, cluster) gives you consistent effect size and prevalence measures to compute 
SpecScore for the already-conserved set. This keeps your specificity ranking aligned to your integration/normalization choices and to the clusters exactly as they were defined.

#### Specificity Score Calculation / Mathematical concept
##### The quantities & formulas
Consider a gene 𝑔 and a cluster 𝑐. Let “others” (not 𝑐) be ¬𝑐.

1) Average expression & log fold-change
   Seurat’s DE returns an effect size as log2 fold change:
   
<img width="189" height="60" alt="Weixin Image_20250909172325_139_103" src="https://github.com/user-attachments/assets/1ef2f6e0-f215-4165-8bb1-b59024445cd8" />


<img width="400" height="70" alt="Weixin Image_20250909172418_140_103" src="https://github.com/user-attachments/assets/acf22261-4846-4625-b709-8917f3aa51de" />


2) Detection rates (prevalence)
   
   <img width="450" height="40" alt="Weixin Image_20250909172510_141_103" src="https://github.com/user-attachments/assets/f7fa7f51-a169-4e73-a79b-0cdd4d4f3e08" />
   
These are fractions of cells with nonzero expression (after normalization on the chosen assay, e.g. RNA or SCT).

3) Specificity score (my algorithm)
   
<img width="562" height="60" alt="Weixin Image_20250909172610_142_103" src="https://github.com/user-attachments/assets/7466f16f-8c8d-467d-89bd-fd7a7b0be4e1" />

<img width="672" height="120" alt="Weixin Image_20250909172642_143_103" src="https://github.com/user-attachments/assets/33b3a21d-a886-4fd8-8a59-264397dfad7c" />

4) Canonical selection rule
   After computing SpecScore, you apply two simple gates before ranking:
   <img width="407" height="40" alt="image" src="https://github.com/user-attachments/assets/0d6557db-91b8-4c65-a98b-c89732c6b3a3" />
   
Then you take the top 4 (for exmaple) by SpecScore per cluster:

<img width="374" height="40" alt="image" src="https://github.com/user-attachments/assets/a0760f60-6b4d-44fd-b2e3-9d249cb6d5fe" />

<img width="605" height="40" alt="image" src="https://github.com/user-attachments/assets/ff572e70-064a-4c6b-8890-b975cedd325d" />




## Installation
To install `scCanonical` from GitHub, use the following commands in R:

```R
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install scCanonical
devtools::install_github("LingzhangMeng/scCanonical")

# Check
library(scCanonical)
```

## Dependencies
The package requires the following R packages, which should be installed prior to using `scCanonical`:
Seurat, dplyr, ggplot2, reshape2 and ggrepel


## Workflow
The following workflow demonstrates how to use `scCanonical` to process scRNA-seq data, integrate datasets, identify canonical markers, and visualize results.

### 1. Load and Prepare Seurat Objects
Read and visualize individual Seurat objects for different conditions (e.g., Control and Wounded).

```R
# 1. Load packages
library(scCanonical)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)


# Load integrated seurat object
Cell.integrated <- readRDS("//path/Seurat.Integration.rds")
```



### 2. Check Metadata and Cluster Distribution
```R
# Safety checks for required metadata columns
if (!"condition" %in% colnames(Cell.integrated@meta.data)) {
  stop("'condition' metadata column is missing in Cell.integrated. ",
       "Please set it before proceeding, e.g.:\n",
       "  Cell.integrated$condition <- Cell.integrated$your_grouping_column")
}
if (!"seurat_clusters" %in% colnames(Cell.integrated@meta.data)) {
  stop("'seurat_clusters' metadata column is missing. Run clustering first.")
}

# Show cluster x condition breakdown
Idents(Cell.integrated) <- "seurat_clusters"
cat("\nClusters: ", paste(levels(Idents(Cell.integrated)), collapse = ", "), "\n")

```

```
output
Clusters:  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 
```


```R
cat("\nCells per cluster:\n")
print(table(Idents(Cell.integrated)))
```


```
output
    0     1     2     3     4     5     6     7     8     9    10    11    12 
16108  8473  5947  5915  4104  3958  3761  2251  2031  1890  1750  1397  1294 
   13 
 1263 
```

```R
cat("\nCells per cluster x condition:\n")
print(table(Cell.integrated$seurat_clusters, Cell.integrated$condition))
```

```
output
     Para_MVIneg Para_MVIpos Tumor_MVIneg Tumor_MVIpos
  0         5200        3527         3280         4101
  1         1800        5222          443         1008
  2           28         159         5637          123
  3           65          47          435         5368
  4         2444         884          297          479
  5          875         896         1163         1024
  6            2           1            1         3757
  7          974         949          128          200
  8            0           0         2026            5
  9          744        1122           18            6
  10           1           2            0         1747
  11         387         415          436          159
  12         114         352          165          663
  13         417         612           59          175
```

```R
# UMAP visualization (plot not shown; will be inserted later)
DimPlot(Cell.integrated, raster = FALSE, pt.size = 0.5,
        label = TRUE, label.size = 6, label.box = FALSE)
```


<img width="672" height="677" alt="image" src="https://github.com/user-attachments/assets/9d09572d-9c33-4779-a2fd-0225568ca37e" />








## Notes
- Ensure that the Seurat object contains `condition` and `seurat_clusters` metadata columns before running marker identification.
- The `RNA` assay is used by default for differential expression analysis, but you can switch to `SCT` if preferred.
- Adjust parameters like `min.pct`, `logfc.threshold`, and `resolution` based on your dataset.

## Contributing
Contributions are welcome! Please submit issues or pull requests to the [GitHub repository](https://github.com/LingzhangMeng/scCanonical).

## License
This package is licensed under the MIT License. See the LICENSE file for details.

## Contact
For questions or support, please contact the package author via GitHub issues.
