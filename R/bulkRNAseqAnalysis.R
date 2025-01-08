#' Perform Bulk RNA-seq Analysis with Gene Mapping
#'
#' @description
#' This function performs comprehensive bulk RNA-seq analysis including differential
#' expression analysis using DESeq2 and gene mapping from Ensembl IDs to gene symbols.
#'
#' @param work_dir Character string specifying working directory path
#' @param metadata_file Character string specifying path to metadata file
#' @param gtf_file Character string specifying path to GTF file
#'
#' @return List containing:
#' \itemize{
#'   \item results - DESeq2 results table
#'   \item deg - Differentially expressed genes
#'   \item summary - Analysis summary statistics
#' }
#'
#' @examples
#' \dontrun{
#' data_dir <- system.file("extdata", package = "GENE-MAPPING")
#' results <- bulkRNAseqAnalysis(
#'   work_dir = data_dir,
#'   metadata_file = file.path(data_dir, "metadata.csv"),
#'   gtf_file = file.path(data_dir, "annotation.gtf.gz")
#' )
#' }
#'
#' @import DESeq2
#' @import rtracklayer
#' @importFrom stats aggregate
#' @importFrom utils read.csv read.table write.csv
#'
#' @export
bulkRNAseqAnalysis <- function(work_dir, metadata_file, gtf_file) {
  # Input validation
  if (!dir.exists(work_dir)) {
    stop("Working directory does not exist: ", work_dir)
  }
  if (!file.exists(metadata_file)) {
    stop("Metadata file does not exist: ", metadata_file)
  }
  if (!file.exists(gtf_file)) {
    stop("GTF file does not exist: ", gtf_file)
  }
  
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
  tryCatch({
    sample_info <- read.csv(metadata_file, row.names=1)
    if (!"Condition" %in% colnames(sample_info)) {
      stop("Metadata file must contain a 'Condition' column")
    }
    print(head(sample_info))
  }, error = function(e) {
    stop("Error reading metadata file: ", e$message)
  })
  
  # Read count files
  tryCatch({
    file_list <- list.files(pattern="*.txt")
    if (length(file_list) == 0) {
      stop("No count files (*.txt) found in working directory")
    }
    count_list <- lapply(file_list, function(x) {
      counts <- as.matrix(read.table(x, row.names=1, sep="\t"))
      if (any(is.na(counts))) {
        warning("NA values found in count file: ", x)
      }
      return(counts)
    })
    all_counts <- do.call(cbind, count_list)
    print(head(all_counts))
  }, error = function(e) {
    stop("Error reading count files: ", e$message)
  })
  
  # Pre-filtering
  keep <- rowSums(all_counts >= 0) >= min(table(sample_info$Condition))
  counts_filtered <- all_counts[keep,]
  print(paste("Dimensions after filtering:", paste(dim(counts_filtered), collapse=" x ")))
  
  # Map genes using GTF file
  tryCatch({
    gene_id_to_symbol <- read_gtf_and_map_genes(gtf_file)
    print(paste("Number of genes mapped:", length(gene_id_to_symbol)))
  }, error = function(e) {
    stop("Error mapping genes from GTF file: ", e$message)
  })
  
  # Map Ensembl IDs to Gene Symbols
  rownames(counts_filtered) <- gene_id_to_symbol[rownames(counts_filtered)]
  
  # Handle NA values
  na_indices <- which(is.na(rownames(counts_filtered)))
  if (length(na_indices) > 0) {
    warning(paste("Found", length(na_indices), "unmapped genes"))
    rownames(counts_filtered)[na_indices] <- paste0("NA_GENE_", seq_along(na_indices))
  }
  
  # Aggregate duplicate genes
  df <- data.frame(gene = rownames(counts_filtered), counts_filtered)
  duplicate_genes <- df$gene[duplicated(df$gene) | duplicated(df$gene, fromLast = TRUE)]
  if (length(duplicate_genes) > 0) {
    warning(paste("Aggregating", length(unique(duplicate_genes)), "duplicate genes"))
    print("Duplicate genes:", unique(duplicate_genes))
  }
  
  df_aggregated <- aggregate(. ~ gene, data=df, sum)
  rownames(df_aggregated) <- df_aggregated$gene
  df_aggregated$gene <- NULL
  counts_filtered <- as.matrix(df_aggregated)
  
  # DESeq2 analysis
  tryCatch({
    colnames(counts_filtered) <- rownames(sample_info)
    dds <- DESeqDataSetFromMatrix(
      countData = counts_filtered,
      colData = sample_info,
      design = ~ Condition
    )
    dds <- DESeq(dds, test="Wald")
    res <- results(dds, contrast=c("Condition","TREATED", "CONTROL"))
    
    # Filter significant genes
    res_adj <- results(dds, contrast=c("Condition", "TREATED", "CONTROL"), 
                      alpha=0.05)
    deg <- subset(res, abs(log2FoldChange) > 1 & padj < 0.05)
    
    print(paste("Number of DEGs:", nrow(deg)))
  }, error = function(e) {
    stop("Error in DESeq2 analysis: ", e$message)
  })
  
  # Create results directory and save results
  save_dir <- file.path(work_dir, "results")
  dir.create(save_dir, showWarnings = FALSE)
  
  tryCatch({
    write.csv(as.data.frame(res), 
              file = file.path(save_dir, "DESeq2_Results.csv"))
    write.csv(as.data.frame(deg), 
              file = file.path(save_dir, "DESeq2_DEGs.csv"))
    
    # Generate basic plots
    png(file.path(save_dir, "MA_plot.png"))
    plotMA(res)
    dev.off()
    
    png(file.path(save_dir, "Volcano_plot.png"))
    plot(res$log2FoldChange, -log10(res$padj),
         main="Volcano Plot",
         xlab="log2 Fold Change",
         ylab="-log10 adjusted p-value")
    dev.off()
  }, error = function(e) {
    warning("Error saving results: ", e$message)
  })
  
  # Return results
  return(list(
    results = res,
    deg = deg,
    summary = list(
      total_genes = nrow(counts_filtered),
      deg_count = nrow(deg),
      save_directory = save_dir,
      unmapped_genes = length(na_indices),
      duplicate_genes = length(unique(duplicate_genes))
    )
  ))
}

#' Read GTF File and Map Genes
#'
#' @param gtf_path Path to GTF file
#' @return Named vector mapping Ensembl IDs to gene symbols
#' @keywords internal
read_gtf_and_map_genes <- function(gtf_path) {
  # Import GTF file
  gtf_data <- import(gtf_path)
  
  # Extract gene annotations
  genes <- subset(gtf_data, type == "gene")
  
  # Create mapping data frame
  gene_info <- data.frame(
    gene_id = sapply(mcols(genes)$gene_id, 
                     function(x) strsplit(x, "\\.")[[1]][1]),
    gene_symbol = mcols(genes)$gene_name,
    stringsAsFactors = FALSE
  )
  
  # Verify mapping quality
  if (any(is.na(gene_info$gene_symbol))) {
    warning("Some Ensembl IDs lack gene symbols")
  }
  
  if (any(duplicated(gene_info$gene_id))) {
    warning("Duplicate Ensembl IDs found")
  }
  
  # Create named vector for lookup
  gene_id_to_symbol <- setNames(gene_info$gene_symbol, gene_info$gene_id)
  return(gene_id_to_symbol)
}
