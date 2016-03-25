# Extract eigenvalues from an ordination analysis
#
# @param x an object returned by an ordination function.
#
# @return A numeric vector containing the eigenvalues.

# generic
eigenvalues <- function(x) {
  UseMethod("eigenvalues")
}


eigenvalues.prcomp <- function(x) { x$sdev^2 }

eigenvalues.PCA <- function(x) { x$eig$eigenvalue }

eigenvalues.rda <- function(x) { as.numeric(x$CA$eig) }

eigenvalues.pca <- function(x) { x$eig }

eigenvalues.pcaRes <- function(x) { as.numeric(x@sDev^2) }


eigenvalues.CA <- function(x) { x$eig$eigenvalue }

eigenvalues.correspondence <- function(x) {
  eig <- x$cor^2
  # MASS allows to extract the last dimension which is meaningless (eigenvalue ~ 0)
  # discard it if needed
  eig <- eig[1:npc(x)]
  return(eig)
}

eigenvalues.ca <- function(x) { x$sv^2 }
