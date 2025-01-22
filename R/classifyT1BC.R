#' Classify T1 bladder cancer samples
#'
#' This samples provides class calls and class probabilities given log-transformed
#' normalized expression data. The expression data must be in the form of a
#' matrix of data.frame, where rows are genes and columns are samples.
#'
#' @param gexp a matrix or data.frame with rows as genes and columns as samples.
#' @param id_type a string, representing the type of gene annotation used in the
#'  expression matrix. One of: "symbol", "ensembl_id" or "entrez_id" for
#'  HGNC symbol, ENSEMBL gene ID, and Entrez Gene ID, respectively.
#' @param data_type a string, representing the type of data to be classified.
#' Can be one of: "rnaseq" or "array". This will change the model used to
#' predict the T1 classes.
#' @param model_type a string, one of "sensitive" for most accurate but less
#' robust model (`enet`) or "robust" for less accurate but more robust model (`knn`).
#' @param min_cor a numeric between 0 and 1. Used only for `data_type = 'array'`.
#' If sample correlation to the centroids of every class fall below this threshold,
#' the sample will be unclassified.
#'
#' @return A data.frame with `ncol(gexp)` rows and 8 columns:
#' * `id` contains the colnames of the original expression matrix
#' * `class` contains the class prediction for the sample
#' If data_type = 'rnaseq':
#' * `class_prob` shows the model probability for the assigned class
#' * `cl_1_prob` to `cl_5_prob` contain the probabilities for each class
#' If data_type = 'array':
#' * `class_cor` shows the correlation of the sample to the centroid of the
#' assigned class.
#' * `cl_1_cor` to `cl_5_cor` contain the correlation of the sample to the
#' centroid of each class
#'
#' @export
#' @importFrom stats predict
#' @importFrom magrittr `%>%`
#' @importFrom dplyr rename_all mutate select case_when
#' @importFrom tidyr everything
#' @importFrom tibble rownames_to_column
classifyT1BC <- function(gexp,
                         id_type = c("symbol", "ensembl_id",
                                      "entrez_id")[1],
                         data_type = c("rnaseq", "array")[1],
                         model_type = c("sensitive", "robust")[1],
                         min_cor = 0.3) {
    #-- Get features from genes4classification
    if(! id_type %in% c("symbol", "ensembl_id", "entrez_id")) {
        stop("`id_type` must be one of 'symbol', 'ensembl_id' or 'entrez_id'.")
    }
    if (data_type == "rnaseq") {
        feats <- genes4classification[,id_type]
    } else if(data_type == "array") {
        feats <- centroids[,id_type]
    } else {
        stop("`data_type` must be one of 'rnaseq' or 'array'.")
    }

    #-- Check input data
    if(!any(class(gexp) %in% c("data.frame", "matrix"))) {
        stop("`gexp` must be a data.frame or matrix")
    } else if (any(class(gexp) %in% c("data.frame"))) {
        gexp <- as.matrix(gexp)
    }
    #-- Address missing
    missing <- feats[! feats %in% rownames(gexp)]
    present_feats <- feats[feats %in% rownames(gexp)]

    if(length(missing) >= 300) {
        stop("We could not find compatible features in this dataset. Are you sure you used the correct `id_type`?")
    }

    if(length(missing) > 0) {
        msg <- paste("There are markers missing from this dataset.
Missing genes:", paste(missing, collapse = ", "))
        warning(msg)
    }

    if(data_type == "rnaseq") {
        #-- Impute missing
        if(length(missing) > 0) {
            val <- mean(gexp[present_feats,], na.rm = TRUE)
            mss_mat <- matrix(val, nrow = length(missing), ncol = ncol(gexp),
                              dimnames = list(missing, colnames(gexp)))
            gexp <- rbind(gexp, mss_mat)
        }
        #-- Convert feats to symbol
        if (id_type != "symbol") {
            gexp2 <- gexp %>%
                as.data.frame() %>%
                rownames_to_column("id") %>%
                right_join(select(genes4classification, "symbol", id = !! id_type), by = "id") %>%
                select(-id) %>%
                column_to_rownames("symbol") %>%
                as.matrix()
        } else {
            gexp2 <- gexp
        }
        data <- gexp2[feats,] %>% t()

        #-- Predict
        if(model_type == "sensitive") {
            pred_cls <- predict(t1BC_model, data)
        } else {
            pred_cls <- predict(t1BC_knn, data)
        }
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
    } else {
        message("-- Using model to predict on array data (nearest shrunken centroid)...")
        #-- Make centroid
        cent_mat <- as.matrix(centroids[,1:5])
        rownames(cent_mat) <- centroids[,id_type]
        colnames(cent_mat) <- as.character(1:5)
        cent_mat <- cent_mat[present_feats,]

        #-- Correlate samples
        scors <- cor(cent_mat, gexp[present_feats,]) %>% t() %>% as.data.frame()

        #-- Make return
        class_res <- data.frame(id = rownames(scors),
                                class = as.character(apply(scors, 1, which.max)),
                                class_cor = apply(scors, 1, max))
        colnames(scors) <- paste0("cl_", colnames(scors), "_cor")
        class_res <- cbind(class_res, scors)
        rownames(class_res) <- NULL

        #-- Filter unclassified
        class_res <- class_res %>%
            mutate(class = ifelse(class_cor < min_cor, "unc", class),
                   class_cor = ifelse(class_cor < min_cor, NA, class_cor)) %>%
            mutate(class = factor(class))

    }

    return(class_res)
}
