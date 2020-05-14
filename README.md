
# classifyT1BC: Classification of T1 Bladder Cancer

## Install

Installing requires R version 3.5 or higher and devtools.
    
    if(!require(devtools)) install.packages("devtools")
    devtools::install_github("csgroen/classifyT1BC")

## Usage
The package consists of the `classifyT1BC` function. It requires as input an RNA-seq gene expression matrix of log-transformed, normalized counts. (e.g. log2FPKM). It supports HGNC Symbol, Ensembl Gene ID and Entrez Gene ID.

    class_results <- classifyT1BC(gexp, type_id = "symbol")
