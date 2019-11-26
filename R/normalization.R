# normalization.R

#' Normalization
#'
#' Performs preprocessing on the raw count dataset.
#'
#' @param x The expression matrix that needs preprocessing or log-normalization.
#' @param percent Genes that are expressed in less than 100*percent\% of the cells are filtered out.
#' @return The preprocessed and log-normalized expression matrix.
#'
#' @examples
#' x = c(1, 2, 3, 4, 5, 6, 7, 8)
#' \dontrun{
#' preprocessing(x)
#' }
preprocessing <- function(x, percent){
  n <- dim(x)[2]
  gene.exprs.count <- rowSums(x != 0)
  x <- x[gene.exprs.count > n * percent, ]
  return(x)
}

#' Normalization
#'
#' Performs log-normalization on the raw count dataset.
#'
#' @param x The expression matrix that needs log-normalization.
#' @return The log-normalized expression matrix.
#'
#' @examples
#' x = c(1, 2, 3, 4, 5, 6, 7, 8)
#' \dontrun{
#' log_normalization(x)
#' }
log_normalization = function(x){
  sf <- colSums(x)/stats::median(colSums(x))
  return(log(sweep(x, 2, sf, '/')+1))
}

# [END]
