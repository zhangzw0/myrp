# cell_clustering.R

#' Cell clustering analysis
#'
#' @param count_data The count expression matrix. The rows correspond to genes and
#' the columns correspond to cells. Can be sparse.
#' @param percent Genes that are expressed in less than 100*percent\% of the cells are filtered out. Default is 0.1.
#' @param ncores Number of cores to use. Default is 4.
#' @param perplexity Perplexity parameter of function Rtsne.
#' @param need.imputation Whether the input matrix needs imputation. Default is FALSE.
#' @param imputed.data Whether the input matrix has been imputed. Default is TRUE.
#' @return The t-SNE visualization and aRI value of the dataset.
#'
#' @author Yiqiu Tan
#'
#' @examples
#' data('ipsc_saver')
#' cell_clustering(ipsc_saver)
#' @export
cell_clustering <- function(count_data, percent=0.1, ncores=4, perplexity = 30,
                            need.imputation=FALSE, imputed.data = TRUE) {
  if (!is.matrix(count_data)) {
    stop("count_data must be a matrix object")
  }
  if (percent > 1 | percent < 0) {
    stop("Percent is not between 0 and 1")
  }
  if (!imputed.data) {
    message("Starting preprocessing ...")
    count_data <- preprocessing(count_data, percent = percent)
    message("Done!")

    if (need.imputation) {
      message("Starting log-normalization ...")
      count_data <- log_normalization(count_data)
      message("Done!")
      message("Starting imputation using SAVER ...")
      count_data <- SAVER::saver(count_data, ncores)$estimate
    }
  }

  message("Starting cell clustering ...")
  kcluster <- length(unique(colnames(count_data)))
  labels <- colnames(count_data)
  # using t-SNE and Kmeans to perform cell clustering
  data <- t(log_normalization(count_data))
  data_out <- Rtsne::Rtsne(data, 2, perplexity = perplexity)
  cl_out <- stats::kmeans(data_out$Y[, 1:2], kcluster)
  aRI <- mclust::adjustedRandIndex(labels, cl_out$cluster)
  message("Done!")

  message("Plot t-SNE visualization of the imputed estimate ...")
  graphics::par(mar = c(0, 0, 0, 0), oma = c(4, 6, 4, 1))
  cols <- c("#e78ac3", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
            "#CC79A7","#8dd3c7", "#ffd92f", "#e78ac3", "#bebada", "#80b1d3", "#fc8d62", "#b3de69")
  graphics::plot(data_out$Y[,1:2], col = scales::alpha(cols[as.numeric(as.factor(labels))], 0.4),
       pch = 19, axes = FALSE, frame.plot = TRUE, ann = FALSE, cex = 0.8)
  graphics::legend("bottomright", legend = round(aRI, 2), bty = "n", cex = 1.4, y.intersp = 0.5)
  message("Done!")
}

# [END]
