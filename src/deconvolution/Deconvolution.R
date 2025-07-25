---
title: "Deconvolution"
author: "Kelvin Mo"
date: "2023-06-30"
output: pdf_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


```{r}
devtools::install_github('satijalab/seurat')
devtools::install_github('satijalab/azimuth')
library(Azimuth)
library(celldex)
library(SingleR)
BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)
install.packages("tidyverse")
library(tidyverse)
```


```{r}
DLBCL.big <- readRDS("../DLBCL_big_final.rds")
```



## Monaco Immune Data Reference = bulk RNA-seq samples of sorted immune cell populations
```{r}
library(Seurat)
# Bulk RNA-seq samples of sorted immune cell populations
MI <- MonacoImmuneData()
df_MI <- readRDS("../DLBCL_big_final.rds")
sce <- as.SingleCellExperiment(df_MI)
results <- SingleR(test=sce, ref=MI, assay.type.test=1, labels = MI$label.main)#label.fine
df_MI$predicted_celltype_MI <- results$labels
table(results$labels)
```
```{r}
# Identify the cells with both "CD4+" and "CD8+" as their cell type
t_cells <- df_MI$predicted_celltype_MI == "CD4+ T cells" | df_MI$predicted_celltype_MI == "CD8+ T cells"
# Update the celltype values to "T cells" for the identified cells
df_MI$predicted_celltype_MI[t_cells] <- "T cells"

DimPlot(df_MI, group.by = "predicted_celltype_MI", reduction = "umap", label = FALSE, label.size=2) + ggtitle("Monaco Immune Data")
```
```{r, fig.width=10}
DimPlot(df_MI, group.by = "predicted_celltype_MI", reduction = "umap", label = FALSE, label.size=2, split.by = "predicted_celltype_MI") + ggtitle("Monaco Immune Data")
```



```{r}
# Define a function to prepare 10X HNSCC atlas data and write a matrix for CIBERSORTx.
write_cibersortx_mx <- function(rna, lm22, metadata, metatype = "cell_type",
                                metasubtype = NULL, codex_full = NULL, split_niche = FALSE) {
  # 1. Prepare granulocyte signatures from LM22 microarray data for CIBERSORTx.
  genes <- sort(intersect(rownames(lm22), rownames(rna)))
  granulocytes <- lm22[genes, c("Neutrophils", "Eosinophils")]
  scale_factor <- mean(colSums(rna[genes, ])) / mean(colSums(granulocytes))
  granulocytes <- granulocytes * scale_factor
  
  # 2. Prepare 10X HNSCC atlas scRNA-seq data for CIBERSORTx.
  cibersortx <- rna[genes, ]
  metadata <- metadata[, c("seurat_clusters", metatype, metasubtype)]
  colnames(metadata)[2] <- "metatype"
  metadata$seurat_clusters <- as.vector(metadata$seurat_clusters)
  metadata$metatype <- as.vector(metadata$metatype)
  
  # 3. Correct cell frequencies by neighborhood cell proportions from CODEX data.
  if (metatype == "neighborhood") {
    if (is.null(metasubtype)) stop("Please provide neighborhood cell subset as cell_type or REMI")
    
    ## Obtain neighborhood cell proportions from CODEX data.
    colnames(metadata)[3] <- "metasubtype"
    metadata$metasubtype <- as.vector(metadata$metasubtype)
    codex_proportions <- codex_full@meta.data[, c("neighborhood", metasubtype, "cell_subset")]
    colnames(codex_proportions)[2] <- "metasubtype"
    granulocyte_ids <- which(codex_proportions$metasubtype == "Granulocytes")
    codex_proportions$metasubtype <- as.vector(codex_proportions$metasubtype)
    codex_proportions$metasubtype[granulocyte_ids] <- codex_proportions$cell_subset[granulocyte_ids]
    codex_proportions <- subset(codex_proportions, metasubtype != "Basophils")
    codex_proportions$metasubtype <- factor(codex_proportions$metasubtype)
    codex_proportions <- codex_proportions %>% group_by(neighborhood, metasubtype, .drop = FALSE) %>% count()
    codex_proportions <- codex_proportions %>% group_by(neighborhood) %>% mutate(total = sum(n))
    codex_proportions <- codex_proportions %>% group_by(neighborhood) %>% mutate(percent = 100 * n / total)
    codex_proportions$n_sample <- round(codex_proportions$percent * 10)
    write.csv(codex_proportions, paste0("CODEX-proportions-by-", metasubtype, ".csv"), row.names = FALSE)
    codex_proportions$neighborhood <- as.vector(codex_proportions$neighborhood)
    codex_proportions$metasubtype <- as.vector(codex_proportions$metasubtype)
    
    ## Create a vector of cell IDs with corrected proportions.
    all_ids <- c()
    for (neighborhood in sort(unique(metadata$metatype))) {
      for (cell_class in unique(metadata$metasubtype)) {
        ids_sample <- which(metadata$metatype == neighborhood & metadata$metasubtype == cell_class)
        if (length(ids_sample) < 10) ids_sample <- which(metadata$metasubtype == cell_class)
        n_sample <- codex_proportions$n_sample[which(codex_proportions$neighborhood == neighborhood & codex_proportions$metasubtype == cell_class)]
        sampled_ids <- sample(ids_sample, size = n_sample, replace = TRUE)
        all_ids <- c(all_ids, sampled_ids)
      }
    }
    
    ## Correct cell frequencies.
    cibersortx <- cibersortx[, all_ids]
    metadata <- metadata[all_ids, ]
    
    ## Add granulocytes.
    #for (neighborhood in sort(unique(metadata$metatype))) {
    #for (granulocyte in c("Neutrophils", "Eosinophils")) {
    #n_sample <- codex_proportions$n_sample[which(codex_proportions$neighborhood == neighborhood & codex_proportions$metasubtype == granulocyte)]
    #metadata <- rbind(metadata, data.frame(seurat_clusters = 20, metatype = neighborhood, metasubtype = rep(granulocyte, n_sample)))
    #cibersortx <- cbind(cibersortx, matrix(data = rep(granulocytes[, granulocyte], n_sample), nrow = nrow(cibersortx), ncol = n_sample))
    #}
    #}
  } else {
    ## Create a vector of sampled cell IDs (max size of all files in CIBERSORTx is 1 GB).
    all_ids <- c()
    for (cell_class in unique(metadata$metatype)) {
      sampled_ids <- which(metadata$metatype == cell_class)
      if (length(sampled_ids) > 1000) sampled_ids <- sample(sampled_ids, size = 1000, replace = FALSE)
      all_ids <- c(all_ids, sampled_ids)
    }
    
    ## Correct cell frequencies.
    cibersortx <- cibersortx[, all_ids]
    metadata <- metadata[all_ids, ]
    
    ## Add granulocytes.
    #for (granulocyte in c("Neutrophils", "Eosinophils")) {
    #  n_sample <- 10  # 10 to avoid filtering out for n = 5 threshold.
    #	 metadata <- rbind(metadata, data.frame(seurat_clusters = 20, metatype = rep(granulocyte, n_sample)))
    #	 cibersortx <- cbind(cibersortx, matrix(data = rep(granulocytes[, granulocyte], n_sample), nrow = nrow(cibersortx), ncol = n_sample))
    #}
  }
  
  ## Write final matrix.
  metadata <- t(metadata)
  if (metatype == "neighborhood") {
    if (split_niche) {
      colnames(cibersortx) <- metadata[3, ]
      # Split T and Malignant cells only.
      if (metasubtype == "REMI") {
        split_ids <- c(which(metadata[3, ] %in% c("Tcell_0", "Tcell_1", "Tcell_2")),
                       which(metadata[3, ] %in% c("Malignant_0", "Malignant_1") & metadata[2, ] %in% c("3", "4", "7", "8")))
      } else if (metasubtype == "cell_type") {
        split_ids <- c(which(metadata[3, ] == "T cells"),
                       which(metadata[3, ] == "Malignant cells" & metadata[2, ] %in% c("3", "4", "7", "8")))
      } else {
        stop("Metasubtype should be REMI or cell_type")
      }
      colnames(cibersortx)[split_ids] <- paste0(metadata[3, split_ids], "_in_", metadata[2, split_ids])
    } else {
      colnames(cibersortx) <- paste0("Niche", metadata[2, ])
    }
  } else {
    colnames(cibersortx) <- metadata[2, ]
  }
  names(dimnames(cibersortx)) <- c("GeneSymbol", "")
  write.table(cibersortx, file = paste0("scRNAseq-for-CIBERSORTx-", metatype,
                                        ifelse(!is.null(metasubtype), paste0("-by-", metasubtype), ""),
                                        ifelse(split_niche, "-split", ""), ".txt"),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}





########################## Process Data ########################################
# Prepare 10X DLBCL atlas scRNA-seqq data with neutrophil gene expression signature
# from LM22 microarray data for CIBERSORTX
DLBCL_atlas <- df_MI
## Import CIBERSORT neutrophil gene expression signature from LM22 matrix.
lm22 <- as.matrix(read.delim("../LM22_source_GEPs.txt", row.names = 1))

## Prepare 10X HNSCC atlas data for CIBERSORTx.
DefaultAssay(DLBCL_atlas) <- "SCT"
rna <- GetAssayData(DLBCL_atlas, assay = "SCT", slot = "data")
rna[which(rna < 0)] <- 0

## Add neutrophil gene expression signature and write data for CIBERSORTx
## separated by cell type or REMI cluster.
write_cibersortx_mx(rna, lm22, metadata = DLBCL_atlas@meta.data, metatype = "predicted_celltype_MI")
```





# CIBERSORTx Workflow
## Tutorial 1: Build a Signature Matrix File from Single-Cell RNA Sequencing Data
## Tutorial 2 - Impute Cell Fractions
## Tutorial 5 - Impute Gene Expression, High-Resolution Mode