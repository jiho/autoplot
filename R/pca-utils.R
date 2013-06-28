# Match the type of element to extract from a PCA object allowing the vocabulary from different packages to co-exist
match.type <- function(type, ...) {

  # define synonyms
  obsTypes <- c("observations", "individuals", "sites", "lines")
  varTypes <- c("variables", "descriptors", "species", "columns")

  # match argument
  type <- match.arg(type, c(obsTypes, varTypes, ...), several.ok=TRUE)

  # reduce synonyms
  type[type %in% obsTypes]  <- "obs"
  type[type %in% varTypes] <- "var"
  type <- unique(type)

  return(type)
}


#' Extract eigenvalues from a Principal Component Analysis object
#'
#' @param x an object resulting from a PCA
#'
#' @param ... pass-through argument
#'
#' @return
#' A numeric vector containing the eigenvalues.
#'
#' @author Jean-Olivier Irisson \email{irisson@@normalesup.org}
#'
#' @examples
#' pca <- prcomp(USArrests, scale = TRUE)
#' eigenvalues(pca)
#'
#' @export
eigenvalues <- function(x, ...) {
  # Generic for the extraction of eigenvalues from a Principal Component Analysis object
  UseMethod("eigenvalues")
}

#' @method eigenvalues prcomp
#' @rdname eigenvalues
#' @export
eigenvalues.prcomp <- function(x, ...) { x$sdev^2 }
#' @method eigenvalues PCA
#' @rdname eigenvalues
#' @export
eigenvalues.PCA <- function(x, ...) { x$eig$eigenvalue }
#' @method eigenvalues pcaRes
#' @rdname eigenvalues
#' @export
eigenvalues.pcaRes <- function(x, ...) { as.numeric(x@sDev^2) }
#' @method eigenvalues pca
#' @rdname eigenvalues
#' @export
eigenvalues.pca <- function(x, ...) { x$eig }
#' @method eigenvalues rda
#' @rdname eigenvalues
#' @export
eigenvalues.rda <- function(x, ...) { as.numeric(x$CA$eig) }


# Generic and methods for the extraction of raw, unscaled, variables scores from a Principal Component Analysis object
var.scores <- function(x, ...) {
  UseMethod("var.scores")
}
var.scores.prcomp <- function(x, ...) { x$rotation }
var.scores.PCA <- function(x, eig, ...) {
  res <- data.frame()
  kind <- c()
  for (cat in c("var", "quanti.sup")) {
    if (!is.null(x[[cat]])) {
      X <- as.data.frame(t(t(x[[cat]]$coord) / sqrt(eig)))
      res <- rbind(res, X)
      kind <- c(kind, rep(cat, nrow(X)))
    }
  }
  attr(res, "kind") <- kind
  return(res)
}
var.scores.pcaRes <- function(x, ...) { x@loadings }
var.scores.pca <- function(x, ...) { x$c1 }
var.scores.rda <- function(x, ...) { x$CA$v }

# Generic and methods for the extraction of raw, unscaled, observations scores from a Principal Component Analysis object
obs.scores <- function(x, ...) {
  UseMethod("obs.scores")
}
obs.scores.prcomp <- function(x, eig, ...) { t(t(x$x) / sqrt(nrow(x$x) - 1) / sqrt(eig)) }
obs.scores.PCA <- function(x, eig, ...) {
  res <- data.frame()
  kind <- c()
  for (cat in c("ind", "ind.sup", "quali.sup")) {
    if (!is.null(x[[cat]])) {
      X <- as.data.frame(t(t(x[[cat]]$coord) / sqrt(eig)) * sqrt(1 / nrow(x[[cat]]$coord)))
      res <- rbind(res, X)
      kind <- c(kind, rep(cat, nrow(X)))
    }
  }
  kind <- gsub("ind", "obs", kind, fixed=TRUE)
  attr(res, "kind") <- kind
  return(res)
}
obs.scores.pcaRes <- function(x, eig, ...) { t(t(x@scores) / sqrt(nrow(x@scores) - 1) / sqrt(eig)) }
obs.scores.pca <- function(x, ...) { x$l1 * sqrt(1 / nrow(x$li)) }
obs.scores.rda <- function(x, ...) { x$CA$u }

# Generic and methods for the extraction of data from a Principal Component Analysis object
get.data <- function(x, ...) {
  # Generic for the extraction of the original data from a Principal Component Analysis object
  UseMethod("get.data")
}
get.data.prcomp <- function(x, ...) { NULL }
get.data.PCA <- function(x, ...) { x$call$X }
get.data.pcaRes <- function(x, ...) { x@completeObs }
get.data.pca <- function(x, ...) { x$tab }
get.data.rda <- function(x, ...) {
  x <- x$CA$Xbar
  center <- attr(x, "scaled:center")
  scale <- attr(x, "scaled:scale")
  unscaledx <- scale(x, center=(-center/scale), scale=1/scale)
  attr(unscaledx, "scaled:center") <- attr(unscaledx, "scaled:scale") <- NULL
  return(unscaledx)
}
