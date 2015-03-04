#' @title Extract information from a Principal Component Analysis
#'
#' Extract information from a Principal Component Analysis such as the scores on the principal components and some fit statistics into a data.frame.
#'
#' @param model an object resulting from a PCA
#'
#' @template pca_params
#'
#' @param ... pass-through argument
#'
#' @return
#' A data.frame containing the original data, when it is supplied or can be extracted from the object, and the additional columns
#'   \item{.id}{the identifier of the observation (row name) or variable (column name)}
#'   \item{.PC#}{the scores of observations or variables on the extracted principal components}
#'   \item{.cos2}{the squared cosine summed over all extracted PCs}
#'   \item{.contrib}{the contribution to the selected PCs}
#'   \item{.kind}{the nature of the data extracted : observations, variables and possibly their status (active or supplementary)}
#'
#' @author Jean-Olivier Irisson \email{irisson@@normalesup.org}
#'
#' @seealso \code{\link{autoplot_pca}} to produce plots based on the output of \code{fortify}.
#' @template pca_seealso
#'
#' @examples
#' # PCA with stats::prcomp
#' pca <- prcomp(USArrests, scale = TRUE)
#' # extract scores of the observations
#' head(fortify(pca))
# # and of the variables
#' fortify(pca, type = "variables")
#' # data is not containe in the `prcomp` object but can be provided
#' head(fortify(pca, data = USArrests, type = "observations"))
#' # select different principal components
#' fortify(pca, type = "var", PC=c(1,3))
#'
#' \dontrun{
#' # PCA with FactoMineR::PCA
#' library("FactoMineR")
#' # add a missing value
#' d <- USArrests
#' d[1,2] <- NA
#' # use supplementary observations and variables
#' pca <- PCA(d, scale = TRUE, graph = FALSE, ind.sup = 2, quanti.sup = 4)
#' # the missing value is replaced by the column mean in the PCA object
#' # the supplementary observation is identified as such
#' head(fortify(pca, data = d))
#' head(fortify(pca))
#' # the supplementary variable is identified as such
#' fortify(pca, type = "variables")
#'
#' # PCA with ad4::dudi.pca
#' library("ade4")
#' pca <- dudi.pca(USArrests, scannf=FALSE)
#' head(fortify(pca))
#' head(fortify(pca, type = "variables"))
#' 
#' # PCA with vegan::rda
#' library("vegan")
#' pca <- rda(USArrests, scale = TRUE)
#' # can use vegan's naming convention
#' head(fortify(pca, type = "sites"))
#' head(fortify(pca, type = "species"))
#'
#' # PCA with pcaMethods::pca, from bioconductor
#' library("pcaMethods")
#' pca <- pca(d, method = "nipals", scale = "uv", completeObs = TRUE, nPcs = 4)
#' # the missing value is imputed by iterative PCA
#' head(fortify(pca))
#' fortify(pca, type = "var")
#'
#' }
#'
#' @importFrom plyr join
#' @name fortify_pca
NULL

fortify_pca <- function(model, data=NULL, type="observations", PC=c(1,2), scaling="auto", ...) {

  # Check arguments
  # type of 
  type <- match.type(type)

  # select the scaling
  if (scaling == "auto") {
    # if automatic, select a scaling appropriate for the type of scores
    scaling <- type
  } else {
    # otherwise select a type of scaling
    if (is.numeric(scaling)) {
      # scaling can be a number, in a way compatible with most of the code out there (vegan in particular)
      if (! scaling %in% 0:3) {
        stop("Scaling should be 0, 1, 2, or 3")
      }
      scaling <- c("none", "obs", "var", "both")[scaling+1]

    } else if (is.character(scaling)) {
      # or a character string, specifying in which "direction" to scale
      scaling <- match.type(scaling, c("raw", "none"))
    
    } else {
      stop("Scaling should be a number between 0 and 3 or a character string")
    }
  }

  # extract eigenvalues
  eig <- eigenvalues(model)
  if (any(PC > length(eig))) {
    stop("At least one of the principal components does not exist")
  }
  if (length(PC) < 1) {
    stop("You must choose at least one principal component")
  }
  # if (length(PC) > 2) {
  #   warning("Extracting information for more than two principal components. The plot might be difficult to read.")
  # }

  # Compute variance explained
  explainedVar <- eig / sum(eig)
  explainedVar <- explainedVar[PC]

  # Extract scores
  if (type == "obs") {
    scores <- obs.scores(model, eig)
    n <- nrow(scores)
  } else {
    # NB: we still need observations scores to compute n, the number of data points
    # TODO this is not very efficient, consider something else
    scores <- obs.scores(model, eig)
    n <- nrow(scores)
    scores <- var.scores(model, eig)
  }
  # extract kind of data (obs/var/supplementary) from the attributes of scores
  kind <- attr(scores, "kind")
  
  # scale scores to compute squared cosine
  scoresScaled <- t(t(scores) * (eig/sum(eig))^(1/2))

  # squared cosine: quality of the representation in the current space
  .cos2 <- ( scoresScaled / sqrt(rowSums(scoresScaled^2)) )^2
  .cos2 <- .cos2[,PC]

  # contribution to the current PCs
  .contrib <- scores^2
  .contrib <- .contrib[,PC]

  if ( length(PC) > 1 ) {
    # the squared cos are additive
    .cos2 <- rowSums(.cos2)
    # contributions are scaled by the percentage of variance explained by each PC so that the contribution displayed is the contribution to the total variance projectable in the current space
    .contrib <- apply(.contrib, 1, function(x,v) {sum(x*v)}, explainedVar)
  }

  # Determine scaling
  if (scaling != "none") {
    
    if ( scaling == type ) {
      # scale observations or variables
      scores <- scoresScaled
    } else if ( scaling == "both") {
      # scale both = scaling 3
      scores <- t(t(scores) * (eig/sum(eig))^(1/4))    
    }

    # scale for biplot anyway (that does not change the geometry of the plot)
    const <- ((n - 1) * sum(eig))^(1/4)
    scores <- scores * const  
  }

  # Turn scores into a data.frame (if they are not already)
  scores <- as.data.frame(scores)

  # Homogenise names
  names(scores) <- paste(".PC", 1:ncol(scores), sep="")

  # extract only the PCs we want
  scores <- scores[PC]

  # Add kind
  if ( ! is.null(kind) ) {
    scores$.kind <- kind
  } else {
    scores$.kind <- type
  }
  # remove contributions of non active elements
  .contrib[scores$.kind != type] <- NA
  
  # Compute identifier row identifier
  .id <- row.names(scores)
  if (is.null(.id)) {
    .id <- 1:nrow(scores)
  }

  res <- data.frame(.id, scores, .cos2, .contrib, stringsAsFactors=FALSE)

  # Add original data, only if we are extracting observations
  if (type == "obs") {
    # fetch data from the model object if necessary
    if (is.null(data)) {
      data <- as.data.frame(get.data(model))
    }
    # if we fetched or already have something, join it with the extracted variables
    if (!is.null(data)) {
      # NB: we have to use join here and not cbind, in case there are supplementary observations
      data$.id <- row.names(data)
      res <- join(data, res, by=".id", type="full")
    }
  }

  # store eigen values and variance explained as attributes
  attr(res, "eig") <- eig[PC]
  attr(res, "explained.variance") <- explainedVar

  return(res)
}

#' @method fortify prcomp
#' @rdname fortify_pca
#' @export
fortify.prcomp <- fortify_pca

#' @method fortify PCA
#' @rdname fortify_pca
#' @export
fortify.PCA <- fortify_pca

#' @method fortify pca
#' @rdname fortify_pca
#' @export
fortify.pca <- fortify_pca

#' @method fortify pcaRes
#' @rdname fortify_pca
#' @export
fortify.pcaRes <- fortify_pca

#' @method fortify rda
#' @rdname fortify_pca
#' @export
fortify.rda <- fortify_pca


