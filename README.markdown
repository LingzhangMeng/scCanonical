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
### pre-processing for integration
```R
# assign "condition" to each group for integration
# this package ONLY recognizes condtion" for down-streaming analysis.
-Control@condition <- "Control"
-Tumor@condition <- "Tumor"
````

```R
# Then prepare list
for (i in 1:length(Cell.list)) {
  Cell.list[[i]] <- SCTransform(Cell.list[[i]], verbose = FALSE)
}

Cell.features <- SelectIntegrationFeatures(object.list = Cell.list, nfeatures = 3000)

Cell.list <- PrepSCTIntegration(object.list = Cell.list, anchor.features = Cell.features, 
                                                         verbose = FALSE)

Cell.anchors <- FindIntegrationAnchors(object.list = Cell.list, dims = 1:30, reduction="rpca", anchor.features = Cell.features, 
                                                                normalization.method = "SCT", verbose = F)
Cell.integrated <- IntegrateData(anchorset = Cell.anchors, normalization.method = "SCT", 
                                 dims = 1:30, new.assay.name = "rpca", k.weight = 50, verbose = F)

Cell.integrated <- ScaleData(Cell.integrated, verbose = T)
```

```R
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


<img width="380" height="335" alt="image" src="https://github.com/user-attachments/assets/9d09572d-9c33-4779-a2fd-0225568ca37e" />



### 3. Find All Markers (FindAllMarkers)
```R
cat("\n--- Running FindAllMarkers ---\n")
DefaultAssay(Cell.integrated) <- "RNA"
Idents(Cell.integrated) <- "seurat_clusters"

all.markers <- FindAllMarkers(Cell.integrated, only.pos = TRUE,
                              min.pct = 0.1, logfc.threshold = 0.25,
                              test.use = "wilcox")
cat("All positive markers: ", nrow(all.markers), " rows\n")
```
```
output
All positive markers:  31151  rows
```

```R
# Remove mitochondrial, ribosomal, and heat-shock genes
rm_bad <- grepl("^MT-|^RPL|^RPS|^HSP|^Mt-|^Rpl|^Rps|^Hsp",
                all.markers$gene, ignore.case = TRUE)
all.markers <- subset(all.markers, !rm_bad)
cat("After filtering MT/RP/HSP: ", nrow(all.markers), " rows\n")
```

```
output
After filtering MT/RP/HSP:  30622  rows
```


### 4. Find Conserved Markers Across Conditions
Using get_conserved_for_all() from scCanonical to find markers that are conserved across the four condition groups (Para_MVIneg, Para_MVIpos, Tumor_MVIneg, Tumor_MVIpos).
```R
cat("\n--- Running get_conserved_for_all ---\n")
cons_out <- get_conserved_for_all(
  Cell.integrated,
  grouping.var        = "condition",
  min.pct             = 0.1,
  logfc.threshold     = 0.25,
  pval_meta_cutoff    = 0.05,
  min.cells.per.group = 3
)

cons.condition <- cons_out$results
skipped.info   <- cons_out$skipped

if (nrow(cons.condition) > 0) {
  cons.condition <- cons.condition %>% dplyr::select(cluster, gene, everything())
}
cat("Conserved markers kept: ", nrow(cons.condition), " rows\n")
```

```
output

```

```R
cat("\nSkipped clusters:\n")
if (nrow(skipped.info) > 0) {
  print(skipped.info)
} else {
  cat("  (none — all clusters passed the conserved test)\n")
}
```

```
output
  cluster        group n_cells
1       6  Para_MVIneg       2
2       6  Para_MVIpos       1
3       6 Tumor_MVIneg       1
4       6 Tumor_MVIpos    3757
5       8 Tumor_MVIneg    2026
6       8 Tumor_MVIpos       5
7      10  Para_MVIneg       1
8      10  Para_MVIpos       2
9      10 Tumor_MVIpos    1747
```
Clusters 6, 8, and 10 were skipped because some condition groups had fewer than 3 cells


### 5. Join with FindAllMarkers and Add Specificity

```R
cat("\n--- Joining conserved markers with FindAllMarkers effect sizes ---\n")

if (nrow(cons.condition) > 0) {
  cons.joined <- cons.condition %>%
    dplyr::left_join(
      all.markers %>% dplyr::select(gene, cluster, avg_log2FC, pct.1, pct.2),
      by = c("gene", "cluster")
    ) %>%
    add_specificity()

  cons.joined <- cons.joined %>%
    dplyr::select(cluster, gene, spec_score, avg_log2FC, pct.1, pct.2, everything())
} else {
  cons.joined <- data.frame()
}
cat("cons.joined rows: ", nrow(cons.joined), "\n")
```

```
output
cons.joined rows:  5804
```

### 6. Select Top Canonical Markers (Conserved Clusters Only)
```R
cat("\n--- Selecting top canonical markers from conserved clusters ---\n")

if (nrow(cons.joined) > 0) {
  cons.joined$cluster <- as.character(cons.joined$cluster)

  canonical <- cons.joined %>%
    dplyr::filter(!is.na(spec_score),
                  avg_log2FC > 0.5,
                  (pct.1 - pct.2) >= 0.20) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = spec_score, n = 4, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(as.numeric(cluster), dplyr::desc(spec_score))
} else {
  canonical <- data.frame()
}
cat("Canonical conserved markers: ", nrow(canonical), " rows\n")
```
```
output
Canonical conserved markers:  44  rows
```

```R
if (nrow(canonical) > 0) {
  cat("\nTop canonical markers (conserved clusters):\n")
  print(as.data.frame(
    canonical %>% dplyr::select(cluster, gene, spec_score, avg_log2FC, pct.1, pct.2)
  ))
}
```

```
output
Top canonical markers (conserved clusters):

cluster	gene	spec_score	avg_log2FC	pct.1	pct.2
0	IGSF6	3.1344433	3.952640	0.883	0.090
0	CSF1R	3.0767013	4.558076	0.720	0.045
0	MS4A7	2.7699865	4.128147	0.752	0.081
0	MSR1	2.6967206	4.824187	0.595	0.036
1	GZMA	2.4130901	3.917354	0.704	0.088
1	TRBC2	2.3664642	3.474984	0.825	0.144
1	CD96	2.3580307	4.093803	0.643	0.067
1	KLRB1	2.3235777	3.765928	0.729	0.112
2	HPGD	1.0871869	3.178909	0.386	0.044
2	MAOB	0.6941880	1.846245	0.542	0.166
2	CLU	0.6841027	1.269207	0.971	0.432
2	PRDX6	0.5967491	1.414097	0.942	0.520
3	TDO2	2.7039071	3.709063	0.814	0.085
3	APOA5	2.2416115	3.170596	0.861	0.154
3	PON3	1.9923961	2.782676	0.923	0.207
3	CP	1.9039460	2.633397	0.978	0.255
4	MZB1	4.5410084	5.096530	0.963	0.072
4	DERL3	3.4554404	5.511069	0.650	0.023
4	JCHAIN	3.1470575	4.463911	0.806	0.101
4	ANKRD36BP2	3.1403390	6.639194	0.479	0.006
5	ADGRL4	5.1232898	6.895410	0.755	0.012
5	RAMP3	5.1227033	6.575999	0.799	0.020
5	FLT1	4.7406876	5.911082	0.831	0.029
5	PTPRB	4.6966070	6.304170	0.764	0.019
7	S100A12	4.0335193	5.211265	0.823	0.049
7	CLEC4E	3.0140721	4.117585	0.809	0.077
7	S100A8	2.7825005	3.710001	0.993	0.243
7	SLC11A1	2.6131570	4.095857	0.707	0.069
9	SFRP5	5.4664269	7.888062	0.698	0.005
9	MMP7	5.1419515	7.201613	0.723	0.009
9	CXCL6	4.5667116	7.237261	0.639	0.008
9	CD24	4.2083650	5.150998	0.870	0.053
11	CLEC9A	5.5861741	6.738449	0.842	0.013
11	XCR1	5.1632043	8.381825	0.619	0.003
11	IDO1	5.0782676	5.823701	0.896	0.024
11	WDFY4	3.4895858	4.351104	0.887	0.085
12	LMOD1	7.0358659	9.125637	0.774	0.003
12	MYH11	6.7935573	8.254626	0.831	0.008
12	MYL9	6.3247504	6.602036	0.995	0.037
12	FRZB	6.2861781	7.799228	0.813	0.007
13	MS4A1	4.2726009	5.696801	0.782	0.032
13	BANK1	3.0616032	5.002620	0.653	0.041
13	LINC00926	2.7726115	6.244620	0.454	0.010
13	CD79A	2.4162529	3.935265	0.682	0.068
```

### 7. Rescue Skipped Clusters
Clusters 6, 8, and 10 were skipped in the conserved analysis because they lacked representation in some condition groups. We rescue them using rescue_skipped_clusters(), which selects top markers from the FindAllMarkers results.
```R
cat("\n--- Rescuing skipped clusters ---\n")

if (nrow(skipped.info) > 0) {
  rescued <- rescue_skipped_clusters(
    all_markers   = all.markers,
    skipped       = skipped.info,
    n             = 4,
    min_log2FC    = 0.5,
    min_delta_pct = 0.20
  )

  cat("\nRescued all_genes rows:  ", nrow(rescued$all_genes), "\n")
  cat("Rescued canonical rows:  ", nrow(rescued$canonical), "\n")

  if (nrow(rescued$canonical) > 0) {
    cat("\nRescued canonical markers (from FindAllMarkers):\n")
    print(as.data.frame(
      rescued$canonical %>%
        dplyr::select(cluster, gene, spec_score, avg_log2FC, pct.1, pct.2)
    ))
  }

  # Combine conserved + rescued
  cons.joined <- dplyr::bind_rows(cons.joined, rescued$all_genes)
  canonical   <- dplyr::bind_rows(canonical, rescued$canonical)

  cat("\nAfter rescue:\n")
  cat("  cons.joined total rows: ", nrow(cons.joined), "\n")
  cat("  canonical total rows:   ", nrow(canonical), "\n")
} else {
  cat("No clusters were skipped — nothing to rescue.\n")
}
```
```
output
Rescued 3 skipped cluster(s): 6, 8, 10 (12 canonical markers selected)

Rescued all_genes rows:   6224 
Rescued canonical rows:   12 

Rescued canonical markers (from FindAllMarkers):
   cluster    gene spec_score avg_log2FC pct.1 pct.2
1       10 PLA2G2A   3.051657   4.063459 0.827 0.076
2       10    LEPR   2.517046   3.462236 0.849 0.122
3       10   CFHR5   2.234530   3.129594 0.871 0.157
4       10     HPD   2.057749   2.918793 0.958 0.253
5        6   NOTUM   6.113809   7.383828 0.837 0.009
6        6    ODAM   5.139074   7.524266 0.687 0.004
7        6  SLC1A2   3.919673   5.497438 0.738 0.025
8        6    TBX3   3.748877   5.448949 0.712 0.024
9        8   MT1CP   6.326658   7.416949 0.862 0.009
10       8   MT1DP   5.584028   6.462995 0.885 0.021
11       8    MT1B   5.477911   5.877587 0.973 0.041
12       8    MT1H   5.145896   5.575185 0.958 0.035

After rescue:
  cons.joined total rows:  12028 
  canonical total rows:    56
```


### 8. Prepare Data for Plotting
```R
cat("\n--- Preparing data for plots ---\n")

# Ensure delta_pct exists everywhere
if (!"delta_pct" %in% colnames(cons.joined)) {
  cons.joined <- cons.joined %>% dplyr::mutate(delta_pct = pct.1 - pct.2)
} else {
  na_idx <- is.na(cons.joined$delta_pct)
  if (any(na_idx)) {
    cons.joined$delta_pct[na_idx] <-
      cons.joined$pct.1[na_idx] - cons.joined$pct.2[na_idx]
  }
}

if (!"delta_pct" %in% colnames(canonical)) {
  canonical <- canonical %>% dplyr::mutate(delta_pct = pct.1 - pct.2)
} else {
  na_idx <- is.na(canonical$delta_pct)
  if (any(na_idx)) {
    canonical$delta_pct[na_idx] <-
      canonical$pct.1[na_idx] - canonical$pct.2[na_idx]
  }
}

# Mark canonical genes in the background data
cons.joined$highlight <- ifelse(cons.joined$gene %in% canonical$gene,
                                "canonical", "other")

# Define cluster order (numeric sort) and factor levels
cluster_order <- sort(unique(as.numeric(as.character(cons.joined$cluster))))
cons.joined$cluster <- factor(cons.joined$cluster, levels = cluster_order)
canonical$cluster   <- factor(canonical$cluster,   levels = cluster_order)

# Color palette (up to 40 clusters)
cb_palette <- c(
  "#ed1299", "#09f9f5", "#246b93", "#cc8e12", "#d561dd", "#c93f00",
  "#ddd53e", "#4aef7b", "#e86502", "#9ed84e", "#AB3282", "#CCC9E6",
  "#8249aa", "#99db27", "#DCC1DD", "#ff523f", "#ce2523", "#f7aa5d",
  "#cebb10", "#03827f", "#931635", "#373bbf", "#a1ce4c", "#ef3bb6",
  "#d66551", "#1a918f", "#ff66fc", "#2927c4", "#7149af", "#57e559",
  "#8e3af4", "#f9a270", "#22547f", "#db5e92", "#edd05e", "#6f25e8",
  "#0dbc21", "#280f7a", "#6373ed", "#5b910f"
)
custom_colors <- cb_palette[seq_along(cluster_order)]
names(custom_colors) <- cluster_order

cat("Clusters in final data:  ",
    paste(levels(cons.joined$cluster), collapse = ", "), "\n")
cat("Clusters with canonical: ",
    paste(sort(unique(as.numeric(as.character(canonical$cluster)))), collapse = ", "), "\n")
```
```
output
Clusters in final data:   0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 
Clusters with canonical:  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
```


### 9. Generate Plots
Two types of plots are generated: a faceted wrapped plot (plot_canonicals) and a single-row inline plot (plot_canonicals_inline). Both highlight the canonical markers.
```R
cat("\n--- Generating plot_canonicals() (faceted wrapped) ---\n")
p1 <- plot_canonicals(cons.joined, canonical, custom_colors)
print(p1)
```
<img width="850" height="600" alt="image" src="https://github.com/user-attachments/assets/6d49c1ee-51d8-4005-8656-a82a166ce4e4" />

```R
cat("\n--- Generating plot_canonicals_inline() (single row) ---\n")
p2 <- plot_canonicals_inline(cons.joined, canonical, custom_colors)
print(p2)
```
<img width="800" height="470" alt="image" src="https://github.com/user-attachments/assets/eec9629b-bffb-4ac9-89bc-90af24b7b7ff" />

















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
