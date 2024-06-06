## Bulk RNA-seq Analysis Gene Mapping Tutorial 

  

### Objective: 

In this tutorial, you will learn how to perform bulk RNA-seq analysis using DESeq2 and local gene annotation mapping. I will guide you through setting up your environment, running the analysis, and interpreting the results. 

  

### Steps: 

  

1. **Clone the Repository:** 

   Open a terminal and run the following command to download the project to your local machine: 

   ```sh 

   git clone https://github.com/mdbabumiamssm/GENE-MAPPING.git 

   ``` 

  

2. **Open RStudio:** 

   Navigate to the `GENE-MAPPING` directory in the terminal and open RStudio: 

   ```sh 

   cd GENE-MAPPING 

   rstudio . 

   ``` 

  

3. **Set the Working Directory:** 

   In RStudio, set the working directory to the `GENE-MAPPING` directory: 

   ```r 

   setwd("/path/to/GENE-MAPPING")  # Replace with your actual path 

   ``` 

  

4. **Install Required Packages:** 

   Install the necessary R packages by running: 

   ```r 

   install.packages(c("DESeq2", "AnnotationDbi", "org.Hs.eg.db", "umap", "ggplot2", "gplots", "RColorBrewer", "VennDiagram", "fgsea", "rtracklayer")) 

   ``` 

  

5. **Load and Source the Script:** 

   Open the `bulkRNAseqAnalysis.R` script in RStudio and source it to load the function: 

   ```r 

   source("bulkRNAseqAnalysis.R") 

   ``` 

  

6. **Run the Analysis:** 

   Call the `bulkRNAseqAnalysis` function with the required parameters. Make sure to replace the paths with your actual paths: 

   ```r 

   bulkRNAseqAnalysis( 

     work_dir = "/path/to/GENE-MAPPING", 

     metadata_file = "group1metadata.csv", 

     gtf_file = "/path/to/Homo_sapiens.GRCh38.112.gtf.gz"  # Ensure you provide the correct path to the GTF file 

   ) 

   ``` 

  

### Results: 

The analysis results will be saved in the `results` directory within your working directory. You can find the differential expression results and significant genes in the output files. 

  

### Conclusion: 

By following this tutorial, you have successfully performed a bulk RNA-seq analysis using DESeq2 and mapped gene annotations locally. This process is crucial for understanding gene expression changes under different conditions. 
