
# classifyT1BC: Classification of T1 Bladder Cancer

## Install

Installing requires R version 3.5 or higher and devtools.
    
    if(!require(devtools)) install.packages("devtools")
    devtools::install_github("csgroen/classifyT1BC")

## Usage
The package consists of the `classifyT1BC` function. It requires as input a log-transformed, normalized RNA-seq gene expression matrix (e.g. log2FPKM). It supports HGNC Symbols, Ensembl Gene IDs and Entrez Gene IDs.

    class_results <- classifyT1BC(gexp, type_id = "symbol")
