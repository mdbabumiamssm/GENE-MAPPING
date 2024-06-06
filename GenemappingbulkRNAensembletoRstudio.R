# Function to perform bulk RNA-seq analysis
bulkRNAseqAnalysis <- function(work_dir, metadata_file, gtf_file) {
  # Set working directory
  setwd(work_dir)
  
  # Load required libraries
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(umap)
  library(ggplot2)
  library(gplots)
  library(RColorBrewer)
  library(VennDiagram)
  library(fgsea)
  library(rtracklayer)
  
  # Read metadata
  sample_info <- read.csv(metadata_file, row.names=1)
  print(head(sample_info))
  
  # Read count files
  file_list <- list.files(pattern="*.txt")
  count_list <- lapply(file_list, function(x) as.matrix(read.table(x, row.names=1, sep="\t")))
  all_counts <- do.call(cbind, count_list)
  print(head(all_counts))
  
  # Pre-filtering
  keep <- rowSums(all_counts >= 0) >= min(table(sample_info$Condition))
  counts_filtered <- all_counts[keep,]
  print(dim(counts_filtered))
  
  gtf_path= "gtf_file"
  # Function to read GTF file and create a mapping of Ensembl IDs to gene symbols
  read_gtf_and_map_genes <- function(gtf_path) {
    # Import GTF file using rtracklayer
    gtf_data <- import(gtf_path)
    # Extract gene annotations and create a data frame
    genes <- subset(gtf_data, type == "gene")
    # Create a data frame with Ensembl gene IDs and corresponding gene symbols
    gene_info <- data.frame(
      gene_id = sapply(mcols(genes)$gene_id, function(x) strsplit(x, "\\..*")[[1]][1]),
      gene_symbol = mcols(genes)$gene_name,
      stringsAsFactors = FALSE
    )
    # Create a named vector for easy lookup
    gene_id_to_symbol <- setNames(gene_info$gene_symbol, gene_info$gene_id)
    return(gene_id_to_symbol)
  }
  
  # Call the function and store the result
  gene_id_to_symbol <- read_gtf_and_map_genes(gtf_file)
  
  # Map Ensembl IDs to Gene Symbols in counts_filtered
  rownames(counts_filtered) <- gene_id_to_symbol[rownames(counts_filtered)]
  
  # Replace NA values with "NA_GENE"
  na_indices <- which(is.na(rownames(counts_filtered)))
  rownames(counts_filtered)[na_indices] <- "NA_GENE"
  
  # Aggregating counts for duplicate gene symbols
  df <- data.frame(gene = rownames(counts_filtered), counts_filtered)
  duplicate_genes <- df$gene[duplicated(df$gene) | duplicated(df$gene, fromLast = TRUE)]
  print("Genes that appear multiple times before aggregation:")
  print(unique(duplicate_genes))
  
  df_aggregated <- aggregate(. ~ gene, data=df, sum)
  print("Genes that were aggregated:")
  print(unique(duplicate_genes))
  
  rownames(df_aggregated) <- df_aggregated$gene
  df_aggregated$gene <- NULL
  counts_filtered <- as.matrix(df_aggregated)
  
  # DESeq2 analysis
  colnames(counts_filtered) <- rownames(sample_info)
  dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = sample_info, design = ~ Condition)
  dds <- DESeq(dds, test="Wald")
  res <- results(dds, contrast=c("Condition","TREATED", "CONTROL"))
  
  # Multiple testing correction is done here with FDR at 5% by default
  res_adj <- results(dds, contrast=c("Condition", "TREATED", "CONTROL"), alpha=0.05)
  deg <- subset(res, abs(log2FoldChange) > 1 & padj < 0.05 )
  
  # Directory to save DEGs and plots
  save_dir <- file.path(work_dir, "results")
  dir.create(save_dir, showWarnings = FALSE)
  
  # Saving and plotting
  write.csv(as.data.frame(res), file = file.path(save_dir, "DESeq2_Results.csv"))
  write.csv(as.data.frame(deg), file = file.path(save_dir, "DESeq2_DEGs.csv"))
}

