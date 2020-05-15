#' Classify T1 bladder cancer samples
#'
#' This samples provides class calls and class probabilities given log-transformed
#' normalized expression data. The expression data must be in the form of a
#' matrix of data.frame, where rows are genes and columns are samples.
#'
#' @param gexp a matrix or data.frame with rows as genes and columns as samples.
#' @param id_type type of gene annation used in the expression matrix. One of:
#' "symbol", "ensembl_id" or "entrez_id" for HGNC symbol, ENSEMBL gene ID, and
#' Entrez Gene ID, respectively.
#'
#' @return A data.frame with `ncol(gexp)` rows and 8 columns:
#' * `id` contains the colnames of the original expression matrix
#' * `class` contains the class prediction for the sample
#' * `class_prob` shows the model probability for the assigned class
#' * `cl_1_prob` to `cl_2_prob` contain the probabilities for each class
#'
#' @export
#' @importFrom stats predict
#' @importFrom magrittr `%>%`
#' @importFrom dplyr rename_all mutate select case_when
#' @importFrom tidyr everything
#' @importFrom tibble rownames_to_column
classifyT1BC <- function(gexp,
                          id_type = c("symbol", "ensembl_id",
                                      "entrez_id")[1]) {
    #-- Get features from genes4classification
    if(! id_type %in% c("symbol", "ensembl_id", "entrez_id")) {
        stop("`id_type` must be one of 'symbol', 'ensembl_id' or 'entrez_id'.")
    }
    feats <- genes4classification[,id_type]

    #-- Check input data
    if(! class(gexp) %in% c("data.frame", "matrix")) {
        stop("`gexp` must be a data.frame or matrix")
    }

    #-- Address missing
    missing <- feats[! feats %in% rownames(gexp)]

    if(length(missing) != 0) {
        if(length(missing) == 300) {
            stop("We could not find compatible features in this dataset. Are you sure you used the correct `id_type`?")
        } else {
            msg <- paste("There are markers missing from this dataset.
Missing genes:", paste(missing, collapse = ", "))
            warning(msg)
        }
        val <- mean(gexp[feats[feats %in% rownames(gexp)],])
        mss_mat <- matrix(val, nrow = length(missing), ncol = ncol(gexp),
                          dimnames = list(missing, colnames(gexp)))
        gexp <- rbind(gexp, mss_mat)
    }
    #-- Rearrange

    #-- Convert feats to symbol
    if (id_type != "symbol") {
        gexp2 <- gexp %>%
            as.data.frame() %>%
            rownames_to_column("id") %>%
            right_join(data, select(genes4classification, "symbol", id = !! id_type), by = "id") %>%
            select(-id) %>%
            column_to_rownames("symbol") %>%
            as.matrix()

    } else {
        gexp2 <- gexp
        }
    data <- gexp[classification_features,] %>% t()


    #-- Make prediction
    pred_cls <- predict(t1BC_model, data)
    probs_cls <- predict(t1BC_model, data, type = "prob") %>%
        rename_all(~ paste0("cl_", ., "_prob"))
    #-- Make return
    class_res <- cbind(class = pred_cls, probs_cls) %>%
        as.data.frame() %>%
        rownames_to_column("id")

    #-- Make assigned class_prob
    class_res <- class_res %>%
        mutate(class_prob = case_when(
            class == 1 ~ cl_1_prob,
            class == 2 ~ cl_2_prob,
            class == 3 ~ cl_3_prob,
            class == 4 ~ cl_4_prob,
            class == 5 ~ cl_5_prob
        )) %>%
        select(id, class, class_prob, everything())
    return(class_res)
}
