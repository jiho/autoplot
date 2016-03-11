#' Tidying methods for a Principal Component Analysis
#'
#' Extract information from a Principal Component Analysis such as the scores on the principal components and some fit statistics.
#'
#' @param x an object returned by a function performing Principal Component Analysis.
#'
#' @param dimensions vector giving the numbers of the principal components to extract. Typically two are extracted to create a plot. 
#'
#' @param data the original dataset, to be concatenated with the output when extracting row scores. When \code{NULL} (the default) data will be extracted from the PCA object when it contains it (i.e. for all cases but \code{\link{prcomp}})
#'
#' @inheritParams scores
#'
#'
#' @param ... pass-through argument
#'
#' @return
#' A data.frame containing the original data, when it is supplied or can be extracted from the object, and the additional columns
#' \describe{
#'   \item{.label}{the identifier of the row or column.}
#'   \item{.PC#}{the scores (coordinates) of the objects on the extracted principal components.}
#'   \item{.cos2}{the squared cosine summed over all extracted PCs, which quantifies the quality of the representation of the object on the extracted PCs.}
#'   \item{.contrib}{the contribution of each object to the selected PCs}
#   \item{.type}{the nature of the data extracted : observations, varia.bles and possibly their status (active or supplementary).}
#' }
#'
# TODO review documentation
# @seealso \code{\link{autoplot_pca}} to produce plots based on the output of \code{augment}.
# @template pca_seealso
#'
# @examples
# # PCA with stats::prcomp
# pca <- prcomp(USArrests, scale = TRUE)
# # extract scores of the observations
# head(augment(pca))
## and of the variables
# augment(pca, type = "variables")
# # data is not containe in the `prcomp` object but can be provided
# head(augment(pca, data = USArrests, type = "observations"))
# # select different principal components
# augment(pca, type = "var", PC=c(1,3))
#
# \dontrun{
# # PCA with FactoMineR::PCA
# library("FactoMineR")
# # add a missing value
# d <- USArrests
# d[1,2] <- NA
# # use supplementary observations and variables
# pca <- PCA(d, scale = TRUE, graph = FALSE, ind.sup = 2, quanti.sup = 4)
# # the missing value is replaced by the column mean in the PCA object
# # the supplementary observation is identified as such
# head(augment(pca, data = d))
# head(augment(pca))
# # the supplementary variable is identified as such
# augment(pca, type = "variables")
#
# # PCA with ad4::dudi.pca
# library("ade4")
# pca <- dudi.pca(USArrests, scannf=FALSE)
# head(augment(pca))
# head(augment(pca, type = "variables"))
# 
# # PCA with vegan::rda
# library("vegan")
# pca <- rda(USArrests, scale = TRUE)
# # can use vegan's naming convention
# head(augment(pca, type = "sites"))
# head(augment(pca, type = "species"))
#
# # PCA with pcaMethods::pca, from bioconductor
# library("pcaMethods")
# pca <- pca(d, method = "nipals", scale = "uv", completeObs = TRUE, nPcs = 4)
# # the missing value is imputed by iterative PCA
# head(augment(pca))
# augment(pca, type = "var")
#
# }
#'
#' @importFrom dplyr select full_join
#' @importFrom broom tidy augment glance
#' @name pca_tidiers
NULL

tidy_pca <- function(x, ...) {
  eig <- eigenvalues(x)
  prop.variance <- eig/sum(eig)
  cum.prop.variance <- cumsum(prop.variance)
  data.frame(
    term=paste0("PC", 1:length(eig)),
    variance=eig,
    prop.variance,
    cum.prop.variance
  )
}

#' @name pca_tidiers
#' @export
tidy.prcomp <- tidy_pca
#' @name pca_tidiers
#' @export
tidy.PCA <- tidy_pca
#' @name pca_tidiers
#' @export
tidy.rda <- tidy_pca
#' @name pca_tidiers
#' @export
tidy.pca <- tidy_pca
#' @name pca_tidiers
#' @export
tidy.pcaRes <- tidy_pca


augment_pca <- function(x, data=NULL, dimensions=c(1,2), type="row", scaling=type, ...) {

  # check arguments
  type <- match_type(type)

  # extract eigenvalues and variance explained
  eig <- tidy(x)
  # and number of dimensions
  n <- nrow(eig)
  if (length(dimensions) < 1) {
    stop("You must choose at least one dimension.")
  }
  if (any(dimensions > n)) {
    stop("At least one of the dimensions requested does not exist.")
  }
  # if (length(dimensions) > 2) {
  #   warning("Extracting information for more than two dimensions. The plot might be difficult to read.")
  # }
  eig <- eig[dimensions,]
  
  # extract scores
  sco <- scores(x, type=type, scaling=scaling)
  sco_num <- sco[,1:n]
  
  # squared cosine: quality of the representation in the current space
  cos2 <- ( sco_num / sqrt(rowSums(sco_num^2)) )^2

  # contribution to each dimension
  contrib <- sco_num^2

  # reduce to the dimensions of interest
  sco <- sco[,c(names(sco)[dimensions], "label")]
  cos2 <- cos2[,dimensions]
  contrib <- contrib[,dimensions]

  # collapse cos2 and contrib in the current space
  if ( length(dimensions) > 1 ) {
    # the squared cos are additive
    cos2 <- rowSums(cos2)
    # contributions are scaled by the proportion of variance explained by each PC so that the contribution displayed is the contribution to the total variance projectable in the current space
    contrib <- apply(contrib, 1, function(x,v) {sum(x*v)}, v=eig$prop.variance)
  }

  # # remove contributions of non active elements
  # contrib[scores$.type != type] <- NA
  
  # prepare result
  sco <- select(sco, label, 1:ncol(sco))
  res <- data.frame(sco, cos2, contrib, stringsAsFactors=FALSE)
  names(res) <- paste0(".", names(res))

  # add original data, only if we are extracting observations
  if (type == "row") {
    if ( is.null(data) ) {
      # fetch data from the object if necessary
      data <- get_data(x)
    } else {
      # otherwise just add labels
      data$.label <- row.names(data)
    }
    # if we fetched or already have something, join it with the extracted variables
    if (!is.null(data)) {
      # NB: we have to use join here and not cbind, in case there are supplementary observations
      res <- full_join(data, res, by=".label")
    }
  }

  # store eigenvalues, variance explained, etc. and  as attributes
  attr(res, "eig") <- eig
  attr(res, "axes.names") <- paste0(eig$term, " (", round(eig$prop.variance * 100), "%)")

  return(res)
}

#' @name pca_tidiers
#' @export
augment.prcomp <- augment_pca

#' @name pca_tidiers
#' @export
augment.PCA <- augment_pca

#' @name pca_tidiers
#' @export
augment.pca <- augment_pca

#' @name pca_tidiers
#' @export
augment.pcaRes <- augment_pca

#' @name pca_tidiers
#' @export
augment.rda <- augment_pca


