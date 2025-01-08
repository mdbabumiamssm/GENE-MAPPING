library(testthat)
library(GENE-MAPPING)

test_that("bulkRNAseqAnalysis input validation works", {
  # Test invalid working directory
  expect_error(
    bulkRNAseqAnalysis(
      work_dir = "nonexistent_dir",
      metadata_file = "metadata.csv",
      gtf_file = "annotation.gtf.gz"
    ),
    "Working directory does not exist"
  )
  
  # Test invalid metadata file
  temp_dir <- tempdir()
  expect_error(
    bulkRNAseqAnalysis(
      work_dir = temp_dir,
      metadata_file = "nonexistent.csv",
      gtf_file = "annotation.gtf.gz"
    ),
    "Metadata file does not exist"
  )
  
  # Test invalid GTF file
  expect_error(
    bulkRNAseqAnalysis(
      work_dir = temp_dir,
      metadata_file = file.path(temp_dir, "metadata.csv"),
      gtf_file = "nonexistent.gtf.gz"
    ),
    "GTF file does not exist"
  )
})

test_that("read_gtf_and_map_genes works correctly", {
  # Create a minimal test GTF file
  test_gtf <- tempfile(fileext = ".gtf")
  writeLines(
    'chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG00000000001"; gene_name "GENE1";',
    test_gtf
  )
  
  # Test function
  result <- read_gtf_and_map_genes(test_gtf)
  
  # Check output structure
  expect_type(result, "character")
  expect_named(result)
  expect_equal(result["ENSG00000000001"], "GENE1")
  
  # Clean up
  unlink(test_gtf)
})
