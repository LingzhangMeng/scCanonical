# scCanonical: An R Package for Canonical Marker Identification in scRNA-seq Analysis

## Overview
The `scCanonical` R package is designed to identify and rank canonical markers for downstream cluster annotations in single-cell RNA sequencing (scRNA-seq) analysis. It is specifically developed to work with integrated Seurat objects created using the `Seurat` R package. The package facilitates the calculation of conserved markers across conditions, computes specificity scores, and visualizes results to aid in the interpretation of scRNA-seq data.

## Installation
To install `scCanonical` from GitHub, use the following commands in R:

```R
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# Install scCanonical
devtools::install_github("LingzhangMeng/scCanonical")
```

## Dependencies
The package requires the following R packages, which should be installed prior to using `scCanonical`:

```R
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggrepel)
```

Load `scCanonical` after installing dependencies:

```R
library(scCanonical)
```

## Workflow
The following workflow demonstrates how to use `scCanonical` to process scRNA-seq data, integrate datasets, identify canonical markers, and visualize results.

### 1. Load and Prepare Seurat Objects
Read and visualize individual Seurat objects for different conditions (e.g., Control and Wounded).

```R
# Load Seurat objects
# Generate seurat objects Control and Wounded in advance with R package Seurat 
Control <- readRDS("Path/Control.rds")
Wounded <- readRDS("Path/Wounded.rds")
```

### 2. Data Integration
Prepare the Seurat objects for integration by assigning condition metadata and performing SCTransform, feature selection, and integration.

```R
# Assign condition metadata
# In this tutorial, we use "condition" as the key word for defining groups 
Control$condition <- "Control"
Wounded$condition <- "Wounded"
```
```R
# These are the standard procedure to create integrated seurat object
# Create a list of Seurat objects
Cell.list <- list(Control, Wounded)

# Remove individual objects to save memory
rm(Control, Wounded, p.Control, p.Wounded)

# Apply SCTransform to each object
for (i in 1:length(Cell.list)) {
  Cell.list[[i]] <- SCTransform(Cell.list[[i]], verbose = FALSE)
}

# Select integration features
Cell.features <- SelectIntegrationFeatures(object.list = Cell.list, nfeatures = 3000)

# Prepare for integration
Cell.list <- PrepSCTIntegration(object.list = Cell.list, anchor.features = Cell.features, verbose = FALSE)

# Find integration anchors
Cell.anchors <- FindIntegrationAnchors(object.list = Cell.list, dims = 1:30, reduction = "rpca", 
                                      anchor.features = Cell.features, normalization.method = "SCT", verbose = FALSE)

# Integrate datasets
Cell.integrated <- IntegrateData(anchorset = Cell.anchors, normalization.method = "SCT", 
                                dims = 1:30, new.assay.name = "rpca", k.weight = 50, verbose = FALSE)

# Clean up
rm(Cell.list, Cell.anchors, Cell.features)

# Perform integrated analysis
Cell.integrated <- ScaleData(Cell.integrated, verbose = TRUE)
Cell.integrated <- RunPCA(Cell.integrated, npcs = 30, verbose = TRUE)
Cell.integrated <- RunUMAP(Cell.integrated, reduction = "pca", dims = 1:30, verbose = TRUE)
Cell.integrated <- FindNeighbors(Cell.integrated, reduction = "pca", dims = 1:30, verbose = TRUE)
Cell.integrated <- FindClusters(Cell.integrated, pc.use = 1:10, resolution = 0.42, group.singletons = TRUE, verbose = FALSE)

# Visualize integrated data
DimPlot(Cell.integrated, raster = FALSE, pt.size = 0.5, label = TRUE, label.size = 6, label.box = FALSE)
```
<img width="889" height="735" alt="Weixin Image_20250909164956_133_103" src="https://github.com/user-attachments/assets/9b55bd17-386e-436e-943b-cc35d161c8f2" />



### 3. Marker Identification
Identify conserved markers across conditions and compute specificity scores.

```R
# Safety checks
if (!"condition" %in% colnames(Cell.integrated@meta.data)) {
  stop("'condition' metadata column is missing in Cell.integrated.")
}
if (!"seurat_clusters" %in% colnames(Cell.integrated@meta.data)) {
  stop("'seurat_clusters' metadata column is missing. Run clustering first.")
}

# Set cluster identities
Idents(Cell.integrated) <- "seurat_clusters"

# Find all positive markers
DefaultAssay(Cell.integrated) <- "RNA"
all.markers <- FindAllMarkers(Cell.integrated, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "wilcox")
cat("All positive markers: ", nrow(all.markers), " rows\n")

# Remove mitochondrial, ribosomal, and heat-shock genes
rm_bad <- grepl("^MT-|^RPL|^RPS|^HSP|^Mt-|^Rpl|^Rps|^Hsp", all.markers$gene, ignore.case = TRUE)
all.markers <- subset(all.markers, !rm_bad)
cat("First-pass markers computed after filtering: ", nrow(all.markers), " rows\n")

# View results
View(all.markers)
```
<img width="665" height="1064" alt="Weixin Image_20250909165253_134_103" src="https://github.com/user-attachments/assets/a1fa092b-1273-41d3-8df6-10a96ddc11ed" />

```R
# Find conserved markers
cons_out <- get_conserved_for_all(
  Cell.integrated,
  grouping.var = "condition",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  pval_meta_cutoff = 0.05,
  min.cells.per.group = 3
)

cons.condition <- cons_out$results
skipped.info <- cons_out$skipped
cons.condition <- cons.condition %>% select(cluster, gene, everything())
cat("Conserved markers kept: ", ifelse(nrow(cons.condition) > 0, nrow(cons.condition), 0), " rows\n")
if (nrow(cons.condition) > 0) View(cons.condition)
if (nrow(skipped.info) > 0) View(skipped.info)
```
<img width="1905" height="860" alt="Weixin Image_20250909165512_135_103" src="https://github.com/user-attachments/assets/b3a3e5e3-ff00-4f97-bc60-2e096b6cada9" />
<img width="270" height="155" alt="Weixin Image_20250909165523_136_103" src="https://github.com/user-attachments/assets/3b198841-fe23-4ea6-ad04-f618ebe192ef" />


```R
# Compute specificity scores
cons.joined <- cons.condition %>%
  dplyr::left_join(
    all.markers %>% dplyr::select(gene, cluster, avg_log2FC, pct.1, pct.2),
    by = c("gene", "cluster")
  ) %>%
  add_specificity()

cons.joined <- cons.joined %>% 
  select(cluster, gene, spec_score, avg_log2FC, pct.1, pct.2, everything())
cat("Conserved + specificity (joined) rows: ", ifelse(nrow(cons.joined) > 0, nrow(cons.joined), 0), "\n")
if (nrow(cons.joined) > 0) View(cons.joined)
```
<img width="2475" height="1075" alt="Weixin Image_20250909165801_137_103" src="https://github.com/user-attachments/assets/8a63ec8b-01d8-411d-9890-3ed92089a1b4" />

```R
# Select top canonical markers
# Here in the below script, you can set either "n = 6" to selct top 6 canonical markers, or set "n = 4" to select top 4 canonical markers
# or to set the other top canonical numbers, for exmaple, n = 1, n = 2, n = 3, n=4, n = 5, n = 6, n = 7.
# I suggest pick a number from 2 to 4, dependes on your purpose
canonical <- cons.joined %>%
  dplyr::filter(!is.na(spec_score), avg_log2FC > 0.5, (pct.1 - pct.2) >= 0.20) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = spec_score, n = 4, with_ties = FALSE) %>%  
  dplyr::ungroup() %>%
  select(cluster, gene, spec_score, avg_log2FC, pct.1, pct.2, everything()) %>%
  mutate(cluster = as.numeric(as.character(cluster))) %>%
  arrange(cluster, desc(spec_score))
cat("Canonical conserved markers kept: ", ifelse(nrow(canonical) > 0, nrow(canonical), 0), " rows\n")
if (nrow(canonical) > 0) View(canonical)
```

<img width="2506" height="266" alt="Weixin Image_20250909170655_138_103" src="https://github.com/user-attachments/assets/75c88eab-3fee-44c6-8bed-480efe39d6d7" />



#### Specificity Score Calculation
The specificity score is computed to rank genes based on their differential expression and specificity to a cluster. The formula is defined as:

\[
\text{spec_score} = \text{avg_log2FC} \times \max(0, \text{pct.1} - \text{pct.2})
\]

Where:
- `avg_log2FC`: Average log2 fold change of the gene in the cluster compared to others.
- `pct.1`: Percentage of cells in the cluster expressing the gene.
- `pct.2`: Percentage of cells outside the cluster expressing the gene.

This formula prioritizes genes with high differential expression and specific expression within the cluster. The function `add_specificity` implements this calculation:

```R
add_specificity <- function(df) {
  if (all(c("avg_log2FC", "pct.1", "pct.2") %in% colnames(df))) {
    df$spec_score <- df$avg_log2FC * pmax(0, df$pct.1 - df$pct.2)
  } else {
    df$spec_score <- NA_real_
  }
  df
}
```

### 4. Visualization
Generate faceted plots to compare conserved and canonical markers.

```R
# Add percentage difference
cons.joined <- cons.joined %>% mutate(delta_pct = pct.1 - pct.2)
canonical <- canonical %>% mutate(delta_pct = pct.1 - pct.2)

# Mark canonical genes
cons.joined$highlight <- ifelse(cons.joined$gene %in% canonical$gene, "canonical", "other")

# Define cluster order and colors
cluster_order <- unique(cons.joined$cluster)
cons.joined$cluster <- factor(cons.joined$cluster, levels = cluster_order)
canonical$cluster <- factor(canonical$cluster, levels = cluster_order)

# Custom color palette
cb_palette <- c("#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00", "#ddd53e", "#4aef7b", 
                "#e86502", "#9ed84e", "#AB3282", "#CCC9E6", "#8249aa", "#99db27", "#DCC1DD", "#ff523f",
                "#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6", 
                "#d66551", "#1a918f", "#ff66fc", "#2927c4", "#7149af", "#57e559" ,"#8e3af4" ,"#f9a270",
                "#22547f", "#db5e92", "#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f")
custom_colors <- cb_palette[1:length(cluster_order)]
names(custom_colors) <- cluster_order

# Generate plots
plot_canonicals(cons.joined, canonical, custom_colors)
plot_canonicals_inline(cons.joined, canonical, custom_colors)
```

[Insert example visualizations for `plot_canonicals` and `plot_canonicals_inline` here]

## Output
- **Data Frames**: The workflow generates data frames (`all.markers`, `cons.condition`, `cons.joined`, `canonical`) that can be viewed using `View()` for manual inspection.
- **Plots**: UMAP visualizations and faceted marker plots are saved as PDF and JPEG files or displayed inline.

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
