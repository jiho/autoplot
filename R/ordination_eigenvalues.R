#' Eigenvalues in an ordination analysis
#' 
#' Extract eigenvalues from an ordination analysis (Principal Components Analysis, Correspondence Analysis, etc.) object in a standard way.
#'
#' @param x an object returned by an ordination function.
#'
#' @return A numeric vector containing the eigenvalues.
#'
#' @template pca_seealso
#' @template ca_seealso
#'
#' @examples
#' # Principal Component Analysis
#' pca <- prcomp(USArrests, scale=TRUE)
#' eigenvalues(pca)
#' if (require("FactoMineR")) {
#'   eigenvalues(PCA(USArrests, graph=F))
#' }
#' if (require("vegan")) {
#'   eigenvalues(rda(USArrests, scale=TRUE))
#' }
#' if (require("ade4")) {
#'   eigenvalues(dudi.pca(USArrests, scannf=FALSE, nf=4))
#' }
#' if (require("pcaMethods")) {
#'   eigenvalues(pca(USArrests, scale="uv", nPcs=4))
#' }
#'
#' # Correspondence analysis
#' clr <- HairEyeColor[,,1]
#' if (require("FactoMineR")) {
#'   eigenvalues(CA(clr, graph=F))
#' }
#' if (require("MASS")) {
#'   eigenvalues(corresp(clr, nf=3))
#' }
#' if (require("ca")) {
#'   eigenvalues(ca(clr))
#' }
#'
#' @export
eigenvalues <- function(x) {
  # Generic for the extraction of eigenvalues from a Principal Component Analysis object
  UseMethod("eigenvalues")
}

#' @name eigenvalues
#' @export
eigenvalues.prcomp <- function(x) { x$sdev^2 }
#' @name eigenvalues
#' @export
eigenvalues.PCA <- function(x) { x$eig$eigenvalue }
#' @name eigenvalues
#' @export
eigenvalues.rda <- function(x) { as.numeric(x$CA$eig) }
#' @name eigenvalues
#' @export
eigenvalues.pca <- function(x) { x$eig }
#' @name eigenvalues
#' @export
eigenvalues.pcaRes <- function(x) { as.numeric(x@sDev^2) }

#' @name eigenvalues
#' @export
eigenvalues.CA <- function(x) { x$eig$eigenvalue }
#' @name eigenvalues
#' @export
eigenvalues.correspondence <- function(x) { x$cor^2 }
#' @name eigenvalues
#' @export
eigenvalues.ca <- function(x) { x$sv^2 }
