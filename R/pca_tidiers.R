#' Tidying methods for a Principal Component Analysis
#'
#' Extract diagnostics, coordinates on the principal components (i.e. rows scores and columns loadings) and some fit statistics from a Principal Component Analysis.
#'
#' @param x an object returned by a function performing Principal Component Analysis.
#'
#' @param data the original dataset, to be concatenated with the output when extracting row scores. When \code{NULL} (the default) data will be extracted from the PCA object when it contains it (i.e. for all functions but \code{\link[stats]{prcomp}}).
#'
#' @param dimensions vector giving the indexes of the principal components to extract. Typically two are extracted to create a plot. 
#'
#' @param which the type of coordinates in the new space to extract: either "rows", "lines", "observations", "objects", "individuals", "sites" (which are all treated as synonyms) or "columns", "variables", "descriptors", "species" (which are, again, synonyms). All can be abbreviated. By default, coordinates of rows are returned. Row coordinates are commonly called 'scores' and column coordinates usually called 'loadings'.
#'
#' @param scaling scaling for the scores. Can be
#' \describe{
#'   \item{"none" (or 0)}{for raw scores,}
#'   \item{"rows" (or 1, or a synonym of "rows")}{to scale row scores by the eigenvalues,}
#'   \item{"columns" (or 2, or a synonym of "columns")}{to scale column scores by the eigenvalues,}
#'   \item{"both" (or 3)}{to scale both row and column scores.}
#' }
#' By default, scaling is adapted to the type of scores extracted (scaling 1 for row scores, scaling 2 for column scores, and scaling 3 when scores are extracted for a biplot).
#' 
#' @details
#' Scaling of scores follows the conventions of package \code{FactoMineR}. In summary, scaling 0 yields unscaled scores, in scaling 1, row scores are multiplied by
#'   \deqn{\sqrt{n \times eig}}{sqrt(n * eig)}
#' where \eqn{n} is the number of active rows in the ordination and \eqn{eig} are the eigenvalues. In scaling 2, column scores are multiplied by
#'   \deqn{\sqrt{eig}}{sqrt(eig)}
#' In scaling 3 both rows and columns are scaled.
#' 
#' @template param_..._ignored
#'
#' @return
#' For \code{tidy}, a data.frame containing the variance (i.e. eigenvalue), the proportion of variance, and the cumulative proportion of variance associated to each principal component.
#' 
#' For \code{augment}, a data.frame containing the original data (when \code{which="rows"} and \code{data} is supplied or can be extracted from the object) and the additional columns:
#' \describe{
#'   \item{.rownames:}{the identifier of the row or column, extracted from the row or column names in the original data.}
#'   \item{.PC#:}{the coordinates of data objects on the extracted principal components.}
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
#' head(augment(pca, which="row"))
#' head(augment(pca, which="col"))
#' # or use your preferred synonym, possibly abbreviated
#' head(augment(pca, which="obs"))
#' head(augment(pca, which="var"))
#' head(augment(pca, which="descriptors"))
#'
#' # data is not contained in the `prcomp` object but can be provided
#' head(augment(pca, data=USArrests, which="row"))
#' # select different principal components
#' augment(pca, which="col", dim=c(2,3))
#'
#' if (require("FactoMineR")) {
#'   pca <- FactoMineR::PCA(USArrests, graph=FALSE, ncp=4)
#'   head(augment(pca, which="individuals"))
#'   head(augment(pca, which="variables"))
#' }
#'
#' if (require("vegan")) {
#'   pca <- vegan::rda(USArrests, scale=TRUE)
#'   # can use vegan's naming convention
#'   head(augment(pca, which="sites"))
#'   head(augment(pca, which="species"))
#' }
#'
#' if (require("ade4")) {
#'   pca <- ade4::dudi.pca(USArrests, scannf=FALSE)
#'   head(augment(pca))
#'   head(augment(pca, which="variables"))
#' }
#'
#' if (require("pcaMethods")) {
#'   pca <- pcaMethods::pca(USArrests, scale="uv")
#'   head(augment(pca))
#'   augment(pca, which="var")
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
#' @import broom
NULL

#' @include ordination_tidy.R
#' @include ordination_data.R
#' @include ordination_scores.R

#' @name pca_tidiers
#' @export
#' @usage tidy(x, ...)
tidy.prcomp <- tidy_ordination
#' @name pca_tidiers
#' @export
#' @usage NULL
tidy.rda <- tidy_ordination
#' @name pca_tidiers
#' @export
#' @usage NULL
tidy.PCA <- tidy_ordination
#' @name pca_tidiers
#' @export
#' @usage NULL
tidy.pca <- tidy_ordination
#' @name pca_tidiers
#' @export
#' @usage NULL
tidy.pcaRes <- tidy_ordination


augment_pca <- function(x, data=NULL, dimensions=c(1,2), which="row", scaling=which, ...) {

  # check arguments
  which <- match_type(which)

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
  
  # compute eigenvalues
  # (used for scaling contribution and for archiving as attribute)
  eig <- tidy(x)
  
  # extract scores
  sco <- scores(x, which=which, scaling=scaling)
  
  # if all potential PCs are kept and scaling is appropriate for the type of scores, compute cos2 and contrib
  if (n == nc(x) & which == scaling) {
    sco_num <- sco[,1:n]
  
    # squared cosine: quality of the representation in the current space
    cos2 <- ( sco_num / sqrt(rowSums(sco_num^2)) )^2

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
      # TODO check this in particular relative to the scaling
    }
  } else {
    cos2 <- NA
    contrib <- NA
  }
  
  # reduce scores to the dimensions of interest
  sco <- sco[,c("rownames", "type", names(sco)[dimensions])]
  eig <- eig[dimensions,]

  # # remove contributions of non active elements
  # contrib[scores$.type != type] <- NA
  
  # prepare result
  res <- data.frame(sco, cos2, contrib, stringsAsFactors=FALSE)
  names(res) <- paste0(".", names(res))

  # add original data, only if we are extracting observations
  if (which == "row") {
    if ( is.null(data) ) {
      # fetch data from the object if necessary
      data <- get_data(x)
    } else {
      # otherwise just add labels
      data$.rownames <- row.names(data)
    }
    # if we fetched or already have something, join it with the extracted variables
    if (!is.null(data)) {
      # NB: we have to use join here and not cbind, in case there are supplementary observations
      res <- dplyr::full_join(data, res, by=".rownames")
    }
  }
  
  # store eigenvalues, variance explained, etc. and  as attributes
  attr(res, "eig") <- eig
  attr(res, "axes.names") <- ordination_axes_titles(eig)
  
  return(res)
}

#' @name pca_tidiers
#' @export
#' @usage augment(x, data=NULL, dimensions=c(1,2), which="row", scaling=which, ...)
augment.prcomp <- augment_pca
#' @name pca_tidiers
#' @export
#' @usage NULL
augment.rda <- augment_pca
#' @name pca_tidiers
#' @export
#' @usage NULL
augment.PCA <- augment_pca
#' @name pca_tidiers
#' @export
#' @usage NULL
augment.pca <- augment_pca
#' @name pca_tidiers
#' @export
#' @usage NULL
augment.pcaRes <- augment_pca


