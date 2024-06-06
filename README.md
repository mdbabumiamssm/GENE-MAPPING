# Bulk RNA-seq Analysis

This repository contains a script for performing bulk RNA-seq analysis using DESeq2 and local gene annotation mapping.

## Prerequisites

Make sure you have the following R packages installed:
- DESeq2
- AnnotationDbi
- org.Hs.eg.db
- umap
- ggplot2
- gplots
- RColorBrewer
- VennDiagram
- fgsea
- rtracklayer

You can install them using:
```r
install.packages(c("DESeq2", "AnnotationDbi", "org.Hs.eg.db", "umap", "ggplot2", "gplots", "RColorBrewer", "VennDiagram", "fgsea", "rtracklayer"))
