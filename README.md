
# classifyT1BC: Classification of T1 Bladder Cancer

## Install

Installing requires R version 3.5 or higher and devtools.
    
    if(!require(devtools)) install.packages("devtools")
    devtools::install_github("csgroen/classifyT1BC")

## Usage
The package consists of the `classifyT1BC` function. It requires as input a log-transformed, FPKM normalized RNA-seq gene expression matrix (i.e. log2(FPKM+1)). It supports HGNC Symbols, Ensembl Gene IDs and Entrez Gene IDs.

    class_results <- classifyT1BC(gexp, id_type = "symbol")

For array data, please set the `data_type` option to `"array"`. This will use an alternative nearest shrunken centroid model, which is more robust to the differences between RNA-seq and array expression spaces. Note that this classifier is less accurate that the RNA-seq classifier (balanced accuracy of ~85%).

    class_results <- classifyT1BC(gexp, id_type = "symbol", data_type = "array")
