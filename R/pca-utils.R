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


#' Match the type of element to extract from a PCA object
#' 
#' Find the type of element we want to extract from a PCA object, allowing the vocabulary from different packages to co-exist
#'
#' @param type the type of element to be extracted: either "observations", "individuals", "sites", "lines" (which are all synonyms) or "variables", "descriptors", "species", "columns" (which are, again, synonyms)
#'
#' @param ... additional types to consider, on top of the usual ones described above
#'
#' @return
#' A vector of characters containing "obs", "var", or one of the user-supplied types
#'
#' @author Jean-Olivier Irisson \email{irisson@@normalesup.org}
#'
#' @examples
#' match.type("observation")
#' match.type("ind")
#' match.type("line")
#' match.type("n", "none")
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


# Generic and methods for the extraction of raw, unscaled, variables scores from a Principal Component Analysis object
var.scores <- function(x, ...) {
  UseMethod("var.scores")
}
var.scores.prcomp <- function(x, ...) { x$rotation }
var.scores.PCA <- function(x, eig, ...) {
  res <- data.frame()
  for (cat in c("var", "quanti.sup")) {
    if (!is.null(x[[cat]])) {
      X <- as.data.frame(t(t(x[[cat]]$coord) / sqrt(eig)))
      X$.kind <- cat
      res <- rbind(res, X)
    }
  }
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
  for (cat in c("ind", "ind.sup", "quali.sup")) {
    if (!is.null(x[[cat]])) {
      X <- as.data.frame(t(t(x[[cat]]$coord) / sqrt(eig)) * sqrt(1 / nrow(x[[cat]]$coord)))
      X$.kind <- cat
      res <- rbind(res, X)
    }
    res$.kind <- gsub("ind", "obs", res$.kind, fixed=TRUE)
  }
  return(res)
}
obs.scores.pcaRes <- function(x, eig, ...) { t(t(x@scores) / sqrt(nrow(x@scores) - 1) / sqrt(eig)) }
obs.scores.pca <- function(x, ...) { x$l1 * sqrt(1 / nrow(x$li)) }
obs.scores.rda <- function(x, ...) { x$CA$u }

#' Extract scores from a Principal Component Analysis object
#' @export
#'
#' @param x an object resulting from a PCA
#'
#' @param scaling whether to scale each projection on a principal component by the square root of the eigenvalues. Scaling observations makes the distance in the plot of observations proportional to the euclidan distance in the original space. Scaling variables relates the angle between variables to the correlations between them (right angle = no correlation), For a biplot, scaling 1 corresponds to scaling observations but not variables; scaling 2 is the reverse
#'
#' @param ... pass-through argument
#'
#' @return
#' A data.frame of scores
#'
#' @author Jean-Olivier Irisson \email{irisson@@normalesup.org}
#'
#' @examples
#' pca <- prcomp(USArrests, scale = TRUE)
#' head(scores(pca))
#' head(scores(pca, type="var"))
#' head(scores(pca, type="obs", scaling=TRUE)
#' head(scores(pca, type="obs", scaling=FALSE)
scores <- function(x, type, scaling=FALSE, ...) {
  type <- match.type(type)
    
  # extract eigenvalues
  eig <- eigenvalues(x)
  
  # extract scores
  if (type == "obs") {
    scores <- obs.scores(x, eig)
    n <- nrow(scores)
  } else {
    # NB: we still need observations scores to compute n, the number of data points
    # TODO this is not very efficient, consider something else
    scores <- obs.scores(x, eig)
    n <- nrow(scores)
    scores <- var.scores(x, eig)
  }
  
  # scale them if needed
  if (scaling) {
    scores <- t(t(scores) * (eig/sum(eig))^(1/2))
  }
  # scale both = scaling 3
  # not implemented yet. Is it really used?
  # scores <- t(t(scores) * (eig/sum(eig))^(1/4)) * const

  # scale for biplot anyway (that does not change the geometry of the plot)
  const <- ((n - 1) * sum(eig))^(1/4)
  scores <- scores * const

  # turn them into a data.frame (if they are not already)
  scores <- as.data.frame(scores)

  # homogenise names
  names(scores) <- paste(".PC", 1:ncol(scores), sep="")

  return(scores)
}

get.data <- function(x, ...) {
  # Generic for the extraction of the original data from a Principal Component Analysis object
  UseMethod("get.data")
}
get.data.prcomp <- function(x, ...) { NULL }
get.data.PCA <- function(x, ...) { x$call$X }
get.data.pcaRes <- function(x, ...) { model@completeObs }
get.data.pca <- function(x, ...) { x$tab }
get.data.rda <- function(x, ...) {
  x <- pcaV$CA$Xbar
  center <- attr(x, "scaled:center")
  scale <- attr(x, "scaled:scale")
  und <- scale(d, center=(-center/scale), scale=1/scale)
  attr(und, "scaled:center") <- attr(und, "scaled:scale") <- NULL
  return(und)
}
