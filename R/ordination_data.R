# Extract data from an ordination analysis
#
# @param x an object returned by an ordination function.
#
# @return A data.frame containing the data or NULL when the data cannot be found.
#
# @template pca_seealso
# @template ca_seealso
#
# @examples
# xS <- prcomp(USArrests, scale=TRUE)
# xV <- vegan::rda(USArrests, scale=TRUE)
# xF <- FactoMineR::PCA(USArrests, graph=F)
# xA <- ade4::dudi.pca(USArrests, scannf=FALSE, nf=4)
# xM <- pcaMethods::pca(USArrests, scale="uv", nPcs=4)
# head(get_data(xS))
# head(get_data(xF))
# head(get_data(xV))
# head(get_data(xA))
# head(get_data(xM))
#
# clr <- HairEyeColor[,,1]
# xF <- FactoMineR::CA(clr, graph=F)
# xM <- MASS::corresp(clr, nf=3)
# xC <- ca::ca(clr)
# head(get_data(xF))
# head(get_data(xM))
# head(get_data(xC))
#
# @export

#' @include ordination_utilities.R

# "user-facing" function that homogenises the output
get_data <- function(x) {
  d <- get_data_(x)
  if ( ! is.null(d) ) {
    # convert to data.frame
    d <- as.data.frame(d)
    # add rownames as a data column and remove them from the d.f
    d$.rownames <- row.names(d)
    row.names(d) <- NULL
  }
  return(d)
}

# internal generic function to actually extract the data
get_data_ <- function(x) { UseMethod("get_data_") }

# methods for PCA
get_data_.prcomp <- function(x) { NULL }
get_data_.PCA <- function(x) { x$call$X }
get_data_.rda <- function(x) { unscale(x$CA$Xbar) }
get_data_.pca <- function(x) { unscale(x$tab, center=x$cent, scale=x$norm) }
get_data_.pcaRes <- function(x) { x@completeObs }

# methods for CA
get_data_.CA <- function(x) { x$call$X }
get_data_.correspondence <- function(x) { x$Freq }
get_data_.ca <- function(x) {
  d <- x$N
  colnames(d) <- x$colnames
  rownames(d) <- x$rownames
  return(d)
}
