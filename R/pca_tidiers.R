#' Tidying methods for a Principal Component Analysis
#'
#' Extract information from a Principal Component Analysis: the scores (i.e. coordinates) on the principal components and some fit statistics.
#'
#' @param x an object returned by a function performing Principal Component Analysis.
#'
#' @param data the original dataset, to be concatenated with the output when extracting row scores. When \code{NULL} (the default) data will be extracted from the PCA object when it contains it (i.e. for all cases but \code{\link{prcomp}}).
#'
#' @param dimensions vector giving the numbers of the principal components to extract. Typically two are extracted to create a plot. 
#'
#' @inheritParams scores
#'
#' @param ... pass-through argument.
#'
#' @return
#' For \code{tidy}, a data.frame containing the variance (i.e. eigenvalue), the proportion of variance, and the cumulative proportion of variance associated to each principal component.
#' 
#' For \code{augment}, a data.frame containing the original data (when \code{type="rows"} and \code{data} is supplied or can be extracted from the object) and the additional columns:
#' \describe{
#'   \item{.label:}{the identifier of the row or column, extracted from the row or column names in the original data.}
#'   \item{.PC#:}{the scores (i.e., coordinates) of data objects on the extracted principal components.}
#'   \item{.cos2:}{the squared cosine, summed over extracted PCs, which quantifies the quality of the representation of each data point in the space of the extracted PCs. NB: \code{cos2} can only be computed when all possible principal components are extracted in the PCA objects; when it is not the case, \code{cos2} is \code{NA}. In several packages, the number of principal components to keep is an argument of the PCA function (and the default is not "all").}
#'   \item{.contrib:}{the contribution of each object to the selected PCs. NB: same comment as for \code{cos2} regarding the number of PCs kept in the PCA object.}
#'   \item{.type:}{the nature of the data extracted : \code{row} or \code{col}.}
#, and possibly their status (active or supplementary).}
#' }
#'
#' @template pca_seealso
#'
#' @examples
#' pca <- prcomp(USArrests, scale = TRUE)
#'
#' tidy(pca)
#'
#' head(augment(pca, type="row"))
#' head(augment(pca, type="col"))
#' # or use your preferred synonym, possibly abbreviated
#' head(augment(pca, type="obs"))
#' head(augment(pca, type="var"))
#' head(augment(pca, type="descriptors"))
#'
#' # data is not contained in the `prcomp` object but can be provided
#' head(augment(pca, data=USArrests, type="row"))
#' # select different principal components
#' augment(pca, type="col", dim=c(1,3))
#'
#' \dontrun{
#' pca <- FactoMineR::PCA(USArrests, graph=FALSE, ncp=4)
#' head(augment(pca, type="individuals"))
#' head(augment(pca, type="variables"))
#'
#' pca <- vegan::rda(USArrests, scale=TRUE)
#' # can use vegan's naming convention
#' head(augment(pca, type="sites"))
#' head(augment(pca, type="species"))
#'
#' pca <- ade4::dudi.pca(USArrests, scannf=FALSE)
#' head(augment(pca))
#' head(augment(pca, type="variables"))
#'
#' pca <- pcaMethods::pca(USArrests, scale="uv")
#' head(augment(pca))
#' augment(pca, type="var")
#' }

# Additional examples
# # PCA with FactoMineR::PCA
# library("FactoMineR")
# # add a missing value
# d <- USArrests
# d[1,2] <- NA
# # use supplementary observations and variables
# pca <- FactoMineR::PCA(d, scale = TRUE, graph = FALSE, ind.sup = 2, quanti.sup = 4)
# # the missing value is replaced by the column mean in the PCA object
#
 # the supplementary observation is identified as such
# # the missing value is imputed by iterative PCA

#'
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
tidy.rda <- tidy_pca
#' @name pca_tidiers
#' @export
tidy.PCA <- tidy_pca
#' @name pca_tidiers
#' @export
tidy.pca <- tidy_pca
#' @name pca_tidiers
#' @export
tidy.pcaRes <- tidy_pca


augment_pca <- function(x, data=NULL, dimensions=c(1,2), type="row", scaling=type, ...) {

  # check arguments
  type <- match_type(type)

  # and number of dimensions
  n <- npc(x)
  if (length(dimensions) < 1) {
    stop("You must choose at least one dimension.")
  }
  if (any(dimensions > n)) {
    stop("At least one of the dimensions requested does not exist.")
  }
  # if (length(dimensions) > 2) {
  #   warning("Extracting information for more than two dimensions. The plot might be difficult to read.")
  # }
  
  # extract scores
  sco <- scores(x, type=type, scaling=scaling)
  
  # if all potential PCs are kept, compute cos2 and contrib
  if (n == nc(x)) {
    sco_num <- sco[,1:n]
  
    # squared cosine: quality of the representation in the current space
    cos2 <- ( sco_num / sqrt(rowSums(sco_num^2)) )^2
    # TODO review computation in particular w/r to scaling

    # contribution to each dimension
    contrib <- sco_num^2
    # TODO review computation in particular w/r to scaling

    # collapse cos2 and contrib in the current space
    cos2 <- cos2[,dimensions]
    contrib <- contrib[,dimensions]
    if ( length(dimensions) > 1 ) {
      # the squared cos are additive
      cos2 <- rowSums(cos2)
      # contributions are scaled by the proportion of variance explained by each PC so that the contribution displayed is the contribution to the total variance projectable in the current space
      contrib <- apply(contrib, 1, function(x,v) {sum(x*v)}, v=eig$prop.variance)
    }
  } else {
    cos2 <- NA
    contrib <- NA
  }
  
  # reduce scores to the dimensions of interest
  sco <- sco[,c("rownames", "type", names(sco)[dimensions])]

  # # remove contributions of non active elements
  # contrib[scores$.type != type] <- NA
  
  # prepare result
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
      res <- dplyr::full_join(data, res, by=".label")
    }
  }

  return(res)
}

#' @name pca_tidiers
#' @export
augment.prcomp <- augment_pca
#' @name pca_tidiers
#' @export
augment.rda <- augment_pca
#' @name pca_tidiers
#' @export
augment.PCA <- augment_pca
#' @name pca_tidiers
#' @export
augment.pca <- augment_pca
#' @name pca_tidiers
#' @export
augment.pcaRes <- augment_pca


