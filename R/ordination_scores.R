# Extract scores from an ordination analysis
#
# @param x an object returned by an ordination function.
#
# @return A data.frame containing the scores, rownames, and score type

#' @include ordination_utilities.R
#' @include ordination_eigenvalues.R
#' @include ordination_dimensions.R

# generic
scores <- function(x, which="rows", scaling=which, ...) {
  UseMethod("scores")
}


# internal functions to actually extract and scale the scores
scores_ <- function(x, which="rows", scaling=which, ...) {
  # check arguments for the type of scores and scaling
  which <- match_type(which)
  scaling <- match_scaling(scaling)
  
  # get eigenvalues and number of active rows to convert scores computed by the various packages to a common "scale 0"
  eig <- eigenvalues(x)
  nr <- nr(x)
  
  # extract and scale scores
  if (which == "row") {
    scores <- row_scores(x, eig, nr)
    if (scaling %in% c("row", "both")) {
      # TODO biplot scaling is likely wrong. But biplots are "wrong" anyway!
      scores <- scale_scores_n(scores, eig, nr)
    }
  } else if (which == "col") {
    scores <- col_scores(x, eig)
    if (scaling %in% c("col", "both")) {
      scores <- scale_scores(scores, eig)
    }
  }

  # convert to a nicely formatted data.frame
  scores <- as.data.frame(scores)
  # homogenise column names
  npc <- npc(x) # NB: PCs are always extracted in order, by all functions
  names(scores)[1:npc] <- paste0("PC", 1:npc)
  # get rownames as a proper data.frame column
  scores$rownames <- row.names(scores)
  row.names(scores) <- NULL
  # set score type
  scores$type <- which
  
  return(scores)
}

# generics to extract scale 0 scores, that dispatch to appropriate methods
col_scores <- function(x, ...) { UseMethod("col_scores") }
row_scores <- function(x, ...) { UseMethod("row_scores") }

# convert to scale 0
unscale_scores <- function(x, eig) { t(t(x) / sqrt(eig)) } 
unscale_scores_n <- function(x, eig, nr) { t(t(x) / sqrt(nr * eig)) }
unscale_scores_n_1 <- function(x, eig, nr) { t(t(x) / sqrt((nr - 1) * eig)) } 

# define methods for row_scores
row_scores.prcomp <- function(x, eig, nr) { unscale_scores_n_1(x$x, eig, nr) }
row_scores.PCA    <- function(x, eig, nr) { unscale_scores_n(x$ind$coord, eig, nr) } # TODO deal with supplementary
row_scores.rda    <- function(x, ...) { x$CA$u }
row_scores.pca    <- function(x, eig, nr) { unscale_scores_n(x$li, eig, nr) }
row_scores.pcaRes <- function(x, eig, nr) { unscale_scores_n_1(x@scores, eig, nr) }

row_scores.CA             <- function(x, eig, ...) { unscale_scores(x$row$coord, eig) }
row_scores.correspondence <- function(x, ...) {
  # MASS allows to extract the last dimension which is meaningless (eigenvalue ~ 0)
  # discard it if needed
  x$rscore[,1:npc(x)]
}
row_scores.ca             <- function(x, ...) { x$rowcoord }

# define methods for col_scores
col_scores.prcomp <- function(x, ...) { x$rotation }
col_scores.PCA    <- function(x, eig) { unscale_scores(x$var$coord, eig) } # TODO deal with supplementary
col_scores.rda    <- function(x, ...) { x$CA$v }
col_scores.pca    <- function(x, ...) { x$c1 }
col_scores.pcaRes <- function(x, ...) { x@loadings }

col_scores.CA             <- function(x, eig) { unscale_scores(x$col$coord, eig) }
col_scores.correspondence <- function(x, ...) {
  # MASS allows to extract the last dimension which is meaningless (eigenvalue ~ 0)
  # discard it if needed
  x$cscore[,1:npc(x)]
}
col_scores.ca             <- function(x, ...) { x$colcoord }


# scale scores
# convention defined somewhat by Numercial Ecology, Legendre & Legendre
#            implemented in FactoMineR
scale_scores <- function(x, eig) { t(t(x) * sqrt(eig)) } 
scale_scores_n <- function(x, eig, nr) { t(t(x) * sqrt(nr * eig)) }
scale_scores_n_1 <- function(x, eig, nr) { t(t(x) * sqrt((nr - 1) * eig)) } 

# methods
scores.prcomp <- scores_
scores.PCA <- scores_
scores.rda <- scores_
scores.pca <- scores_
scores.pcaRes <- scores_

scores.CA <- scores_
scores.correspondence <- scores_
scores.ca <- scores_

