#' Extract scores from an ordination object
#'
#' @inheritParams eigenvalues
#'
#' @param type the scores (i.e. coordinates in the new space) to extract: either "rows", "lines", "observations", "objects", "individuals", "sites" (which are all treated as synonyms) or "columns", "variables", "descriptors", "species" (which are, again, synonyms). All can be abbreviated. By default, scores of rows are returned.
#'
#' @param scaling scaling for the scores. Can be
#' \describe{
#'   \item{"none" (or 0)}{for raw scores,}
#'   \item{"rows" (or 1, or a synonym of "rows")}{to scale row scores by the eigenvalues,}
#'   \item{"column" (or 2, or a synonym of "column")}{to scale column scores by the eigenvalues,}
#'   \item{"both" (or 3)}{to scale both row and column scores.}
#' }
#' By default, scaling is adapted to the type of scores extracted (scaling 1 for row scores and scaling 2 for column scores).
#' 
#' @details
#' Scaling of scores follows the conventions of package \code{vegan}. In summary, except for scaling 0, scores are all multiplied by a constant:
#'	\deqn{c = \sqrt[4]{(n-1) \times \sum{eig}}}{c = sqrt(sqrt((n-1) * sum(eig)))}
#' where \eqn{n} is the number of active rows in the ordination and \eqn{eig} are the eigenvalues (see function \code{\link{eigenvalues}}).
#' In addition, for scaling 1 (resp. 2), row scores (resp. column scores) are multiplied by:
#' \deqn{\sqrt{\frac{eig}{\sum{eig}}}}{sqrt(eig/sum(eig))}
#' For scaling 3, both row and column scores are multiplied by:
#' \deqn{\sqrt[4]{\frac{eig}{\sum{eig}}}}{sqrt(sqrt(eig/sum(eig)))}
#' For details and justification, see \code{vignette("decision-vegan")}.
#' 
# TODO document return and seealso
#' @examples
#' # Principal Component Analysis
#' pca <- prcomp(USArrests, scale=TRUE)
#' head(scores(pca))
#' head(scores(pca, type="columns"))
#' head(scores(pca, type="columns", scaling=0))
#' head(scores(pca, type="columns", scaling=3))
#'
#' sc <- scores(pca)
#' plot(sc$PC1, sc$PC2, asp=1)
#' text(sc$PC1, sc$PC2, label=sc$label, adj=c(-0.1,0.5))
#'
#' if (require("FactoMineR")) {
#'   head(scores(PCA(USArrests, graph=F)))
#' }
#' if (require("vegan")) {
#'   head(scores(rda(USArrests, scale=TRUE)))
#' }
#' if (require("ade4")) {
#'   head(scores(dudi.pca(USArrests, scannf=FALSE, nf=4)))
#' }
#' if (require("pcaMethods")) {
#'   head(scores(pca(USArrests, scale="uv", nPcs=4)))
#' }
#'
#' # Correspondence analysis
#' clr <- HairEyeColor[,,1]
#' if (require("FactoMineR")) {
#'   CA.res <- CA(clr, graph=F)
#'   scores(CA.res, type="row", scaling="none")
#'   scores(CA.res, type="col", scaling="none")
#'   # for a correspondence analysis, scaling="both" makes most sense
#'   scores(CA.res, type="row", scaling="both")
#'   scores(CA.res, type="col", scaling="both")
#' }
#' if (require("MASS")) {
#'   corresp.res <- corresp(clr, nf=3)
#'   scores(corresp.res, type="row", scaling="both")
#'   scores(corresp.res, type="col", scaling="both")
#' }
#' if (require("ca")) {
#'   ca.res <- ca(clr, nf=3)
#'   scores(corresp.res, type="row", scaling="both")
#'   scores(corresp.res, type="col", scaling="both")
#' }
#'
#' @export
scores <- function(x, type="rows", scaling=type, ...) {
  UseMethod("scores")
}

scores_ <- function(x, type="rows", scaling=type, ...) {
  # TODO add selection of axes here?
  # check arguments for the type of scores and scaling
  type <- match_type(type)
  scaling <- match_scaling(scaling, type)
  
  # compute eigenvalues and number of active rows to:
  # - convert scores computed by the various packages to a common "scale 0"
  # - scale scores afterwards
  eig <- eigenvalues(x)
  nr <- nr(x)
  
  # get and scale scores
  if (type == "row") {
    scores <- row_scores(x, eig, nr)
    scaled <- scale_row_scores(scores, eig, nr, scaling)
  } else if (type == "col") {
    scores <- col_scores(x, eig)
    scaled <- scale_col_scores(scores, eig, nr, scaling)
  }

  # convert to a nicely formatted data.frame
  scaled <- as.data.frame(scaled)
  # homogenise column names
  nc <- length(eig)
  names(scaled)[1:nc] <- paste0("PC", 1:nc)
  # get labels as a proper data.frame column
  scaled$label <- row.names(scaled)
  row.names(scaled) <- NULL
  
  return(scaled)
}

#' @name scores
#' @export
scores.prcomp <- scores_
#' @name scores
#' @export
scores.PCA <- scores_
#' @name scores
#' @export
scores.rda <- scores_
#' @name scores
#' @export
scores.pca <- scores_
#' @name scores
#' @export
scores.pcaRes <- scores_
#' @name scores
#' @export
scores.CA <- scores_
#' @name scores
#' @export
scores.correspondence <- scores_
#' @name scores
#' @export
scores.ca <- scores_


## Utility functions for scores extraction and scaling ----


# Determine the type of scores to extract from an ordination object
# This defines synonyms to allow the vocabulary from different packages to co-exist
match_type <- function(type, ...) {
  # define synonyms
  row_types <- c("rows", "lines", "observations", "objects", "individuals", "sites")
  col_types <- c("columns", "variables", "descriptors", "species")

  # allow abbreviation
  type <- match.arg(type, c(row_types, col_types, ...), several.ok=FALSE)

  # reduce synonyms
  type[type %in% row_types] <- "row"
  type[type %in% col_types] <- "col"

  return(type)
}

# Determine the type scaling for scores extracted from an ordination object
# This uses the same synonyms
match_scaling <- function(scaling, type) {
  # otherwise select a type of scaling
  if (is.numeric(scaling)) {
    # scaling can be a number, in a way compatible with most of the code out there (vegan in particular)
    if (! scaling %in% 0:3) {
      stop("Scaling should be 0, 1, 2, or 3")
    }
    scaling <- c("none", "row", "col", "both")[scaling+1]

  } else if (is.character(scaling)) {
    # or a character string, specifying in which "direction" to scale
    scaling <- match_type(scaling, c("none", "both"))
  
  } else {
    stop("Scaling should be a number between 0 and 3 or a character string")
  }
  
  return(scaling)
}


# Get number of active observations (data rows)
# generic
nr <- function(x) { UseMethod("nr") }
# define methods
nr.prcomp <- function(x) { nrow(x$x) }
nr.PCA    <- function(x) { nrow(x$ind$coord) }
nr.rda    <- function(x) { nrow(x$CA$u) }
nr.pca    <- function(x) { nrow(x$li) }
nr.pcaRes <- function(x) { nrow(x@scores) }

nr.CA <- function(x) { nrow(x$row$coord) }
nr.correspondence <- function(x) { nrow(x$rscore) }
nr.ca <- function(x) { nrow(x$rowcoord) } # TODO handle supplementary



# Extract scale 0 scores

# generics to dispatch to methods
col_scores <- function(x, ...) { UseMethod("col_scores") }
row_scores <- function(x, ...) { UseMethod("row_scores") }

# convert to scale 0
unscale_scores <- function(x, eig) { t(t(x) / sqrt(eig)) } 
unscale_scores_n <- function(x, eig, nr) { t(t(x) / sqrt(nr * eig)) }
unscale_scores_n_1 <- function(x, eig, nr) { t(t(x) / sqrt((nr - 1) * eig)) } 

# define methods
row_scores.prcomp <- function(x, eig, nr) { unscale_scores_n_1(x$x, eig, nr) }
row_scores.PCA <- function(x, eig, nr) { unscale_scores_n(x$ind$coord, eig, nr) } # TODO deal with supplementary
row_scores.rda <- function(x, ...) { x$CA$u }
row_scores.pca <- function(x, eig, nr) { unscale_scores_n(x$li, eig, nr) }
row_scores.pcaRes <- function(x, eig, nr) { unscale_scores_n_1(x@scores, eig, nr) }

row_scores.CA <- function(x, eig, ...) { unscale_scores(x$row$coord, eig) }
row_scores.correspondence <- function(x, ...) { x$rscore }
row_scores.ca <- function(x, ...) { x$rowcoord }

# convert to scale 0

# define methods
col_scores.prcomp <- function(x, ...) { x$rotation }
col_scores.PCA <- function(x, eig) { unscale_scores(x$var$coord, eig) } # TODO deal with supplementary
col_scores.rda <- function(x, ...) { x$CA$v }
col_scores.pca <- function(x, ...) { x$c1 }
col_scores.pcaRes <- function(x, ...) { x@loadings }

col_scores.CA <- function(x, eig) { unscale_scores(x$col$coord, eig) }
col_scores.correspondence <- function(x, ...) { x$cscore }
col_scores.ca <- function(x, ...) { x$colcoord }


# Scale scores
# Conventions are those of vegan. See vignette("decision-vegan") for details.
scale_row_scores <- function(x, eig, nr, scaling="none", ...) {
  # define scaling factors
  prop_eig_2 <- sqrt(eig/sum(eig))
  prop_eig_4 <- sqrt(prop_eig_2)
  const <- ((nr-1)*sum(eig))^(1/4)
  
  # perform scaling
  switch(scaling,
    none = x,
    row = t(t(x) * (prop_eig_2 * const)),
    col = x * const,
    both = t(t(x) * (prop_eig_4 * const))
  )
}

scale_col_scores <- function(x, eig, nr, scaling="none", ...) {
  # define scaling factors
  prop_eig_2 <- sqrt(eig/sum(eig))
  prop_eig_4 <- sqrt(prop_eig_2)
  const <- ((nr-1)*sum(eig))^(1/4)
  
  # perform scaling
  switch(scaling,
    none = x,
    row = x * const,
    col = t(t(x) * (prop_eig_2 * const)),
    both = t(t(x) * (prop_eig_4 * const))
  )
}

