# GENE-MAPPING

## Overview

GENE-MAPPING is a Bioconductor package that provides tools for bulk RNA-seq analysis with automatic gene mapping from Ensembl IDs to gene symbols. The package integrates with DESeq2 for differential expression analysis and provides comprehensive visualization tools.

## Features

- Bulk RNA-seq data analysis using DESeq2
- Automatic gene mapping from Ensembl IDs to gene symbols
- Local gene annotation mapping using GTF files
- Comprehensive visualization tools
- Statistical analysis and reporting

## Installation

You can install the package from Bioconductor:

```r
if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("GENE-MAPPING")
```

## Usage

Basic usage example:

```r
library(GENE-MAPPING)

results <- bulkRNAseqAnalysis(
    work_dir = "path/to/working/directory",
    metadata_file = "metadata.csv",
    gtf_file = "annotation.gtf.gz"
)
```

## Documentation

For detailed documentation, please refer to the package vignette:

```r
browseVignettes("GENE-MAPPING")
```

## Contributing

Please feel free to submit issues and pull requests to our GitHub repository.

## License

This package is licensed under the Artistic-2.0 License.

## Author

MD Babu Mia  
Contact: mdbabumia777@gmail.com
