# Generic and methods for the extraction of data from a Principal Component Analysis object
get_data <- function(x, ...) {
  # Generic for the extraction of the original data from a Principal Component Analysis object
  UseMethod("get_data")
}
get_data.prcomp <- function(x, ...) { NULL }
get_data.PCA <- function(x, ...) { x$call$X }
get_data.pcaRes <- function(x, ...) { x@completeObs }
get_data.pca <- function(x, ...) { x$tab }
get_data.rda <- function(x, ...) {
  x <- x$CA$Xbar
  center <- attr(x, "scaled:center")
  scale <- attr(x, "scaled:scale")
  unscaledx <- scale(x, center=(-center/scale), scale=1/scale)
  attr(unscaledx, "scaled:center") <- attr(unscaledx, "scaled:scale") <- NULL
  return(unscaledx)
}
