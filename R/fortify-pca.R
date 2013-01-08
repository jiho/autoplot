#' Extract information from a Principal Component Analysis into a data.frame
#'
#' @param model an object resulting from a PCA, using \code{\link{prcomp}} in package \code{stats}, \code{\link[FactoMineR:PCA]{PCA}} in package \code{factoMineR} or \code{\link[pcaMethods:pca]{pca}} in package \code{pcaMethods} from bioconductor
#'
#' @param data the original data used to compute the PCA, to be concatenated to the output when extracting observations; when \code{NULL}, the default, the data will be extracted from the PCA object when possible (not for \code{\link{prcomp}})
#'
#' @param type the type of data to extract : observations (i.e. rows, individuals) or variables (i.e. columns, descriptors); can be abbreviated
#'
#' @param PC the principal components to extract; two is usual for plotting
#'
#' @param ... pass-through argument
#'
#' @return
#' A data.frame containing the original data when pertinent and possible and the additional elements
#'   \item{.id}{the identifier of the observation (row name) or variable (column name)}
#'   \item{.PC#}{the scores (for observations) or loadings (for variables) on the extracted principal components}
#'   \item{.cos2}{the squared cosine summed over all extracted PCs}
#'   \item{.contrib}{the contribution to the selected PCs}
#'   \item{.kind}{the nature of the data extracted : observations, variables and possibly their status (active or supplementary)}
#'
#' @author Jean-Olivier Irisson \email{irisson@@normalesup.org}
#'
#' @seealso , \code{\link{prcomp}} in package \code{stats}, \code{\link[FactoMineR:PCA]{PCA}} in package \code{factoMineR} or \code{\link[pcaMethods:pca]{pca}} in package \code{pcaMethods} from bioconductor, \link{autoplot_pca} to produce plots based on the output for \code{fortify}
#'
#' @examples
#' # PCA with stats::prcomp
#' pca <- prcomp(USArrests, scale = TRUE)
#'
#' head(fortify(pca))
#' fortify(pca, type = "variables")
#'
#' # data is not containe in the `prcomp` object but can be provided
#' head(fortify(pca, data = USArrests, type = "observations"))
#'
#' # select different principal components
#' fortify(pca, type = "var", PC=c(1,3))
#'
#' \dontrun{
#' # PCA with FactoMineR::PCA
#' library("FactoMineR")
#'
#' # add a missing value
#' d <- USArrests
#' d[1,2] <- NA
#' # use supplementary observations and variables
#' pca <- PCA(d, scale = TRUE, graph = FALSE, ind.sup = 2, quanti.sup = 4)
#'
#' # the missing value is replaced by the column mean
#' # the supplementary observation is identified as such
#' head(fortify(pca))
#' head(fortify(pca, data = d))
#'
#' # the supplementary variable is identified as such
#' fortify(pca, type = "variables")
#'
#' # PCA with pcaMethods::pca, from bioconductor
#' library("pcaMethods")
#' pca <- pca(d, method = "nipals", scale = "uv", completeObs = TRUE, nPcs = 4)
#'
#' # the missing value is imputed by iterative PCA
#' head(fortify(pca))
#'
#' fortify(pca, type="var")
#'
#' # PCA with vegan::rda
#' library("vegan")
#' pca <- rda(USArrests, scale = TRUE)
#'
#' head(fortify(pca))
#' fortify(pca, type="var")
#'
#' }
#'
#' @importFrom plyr join
#' @name fortify_pca
NULL

#' @method fortify prcomp
#' @rdname fortify_pca
#' @export
fortify.prcomp <- function(model, data=NULL, type=c("observations", "variables"), PC=c(1,2), ...) {
  #
  # Method for stats::prcomp
  #
  # NB: those objects do not contain the original data, so data is NULL by default
  #

  # Checks
  if (any(PC > ncol(model$x))) {
    stop("At least one of the principal components does not exist")
  }
  if (length(PC) < 1) {
    stop("You must choose at least one principal component")
  }
  if (length(PC) > 2) {
    warning("Extracting information for more than two principal components. The plot might be difficult to read.")
  }
  type <- match.arg(type)


  # Compute variance explained
  # eigenvalues (scaled to be homogeneous with other packages)
  eig <- model$sdev^2
  # variance explained by each PC
  explainedVar <- eig / sum(eig)
  explainedVar <- explainedVar[PC]


  # Extract appropriate data
  if (type == "variables") {

    # Variables
    # variable identifier
    .id <- row.names(model$rotation)

    # loadings on the PCs
    loadings <- as.data.frame((t(t(model$rotation)*model$sdev)))
    # TODO select PCs here rather than later and save a bit of memory/cpu?
    # NB: loadings are scaled by std deviation, to be consistent with FactoMineR and ADE4, but this will result in a plot different from the one produced by stat::biplot, where loadings are not scaled
    names(loadings) <- paste(".", names(loadings), sep="")

    # squared cosine: quality of the representation on the current space
    .cos2 <- ( loadings / sqrt(rowSums(loadings^2)) )^2
    .cos2 <- .cos2[,PC]

    # contribution to the current PCs
    .contrib <- as.data.frame(t(t(loadings^2) / eig)) * 100
    .contrib <- .contrib[,PC]

    if (length(PC > 1)) {
      # the squared cos are additive
      .cos2 <- rowSums(.cos2)
      # contributions are scaled by the percentage of variance explained by each PC so that the contribution displayed is the contribution to the total variance projectable in the current space
      .contrib <- apply(.contrib, 1, function(x,v) {sum(x*v)}, explainedVar)
    }

    res <- data.frame(.id, loadings[,PC], .cos2, .contrib, .kind=type)

  }
  else if (type == "observations") {

    # Observations
    # observation identifier
    .id <- row.names(model$x)
    if (is.null(.id)) {
      .id <- 1:nrow(model$x)
    }

    # scores (i.e. coordinates) on the PCs
    scores <- as.data.frame(model$x[,PC])
    names(scores) <- paste(".", names(scores), sep="")
    # TODO check why these are different from ind$coord of FactoMineR (the difference is a scale factor which is always the same: 1.11803398874989)

    # square cosine
    .cos2 <- ( model$x^2 / rowSums(model$x^2) )
    .cos2 <- .cos2[,PC]

    # contributions
    .contrib <- as.data.frame( t( t(model$x^2) * (1/nrow(model$x)) / eig ) ) * 100
    .contrib <- .contrib[,PC]

    if (length(PC > 1)) {
      .cos2 <- rowSums(.cos2)
      .contrib <- apply(.contrib, 1, function(x,v) {sum(x*v)}, explainedVar)
    }

    res <- data.frame(.id, scores, .cos2, .contrib, .kind=type, stringsAsFactors=FALSE)
    if (!is.null(data)) {
      data$.id <- row.names(data)
      res <- join(data, res, by=".id", type="full")
    }
  }


  # store variance explained as attribute
  attr(res, "explained.variance") <- explainedVar

  return(res)
}

#' @method fortify PCA
#' @rdname fortify_pca
#' @export
fortify.PCA <- function(model, data=NULL, type=c("observations", "variables"), PC=c(1,2), ...) {
  #
  # Method for FactoMineR::PCA
  #

  # Checks
  if (any(PC > ncol(model$x))) {
    stop("At least one of the principal components does not exist")
  }
  if (length(PC) < 1) {
    stop("You must choose at least one principal component")
  }
  if (length(PC) > 2) {
    warning("Extracting information for more than two principal components. The plot might be difficult to read.")
  }
  type <- match.arg(type)


  # Compute variance explained
  explainedVar <- model$eig$eigenvalue / sum(model$eig$eigenvalue)
  explainedVar <- explainedVar[PC]


  # Extract appropriate data
  # determine what to extract
  if (type == "variables") {
    extract <- c("var", "quanti.sup")
    } else if (type == "observations") {
    extract <- c("ind", "ind.sup", "quali.sup")
  }

  # prepare storage
  d <- list()

  for (i in extract) {

    if (!is.null(model[[i]])) {

      # variable identifier
      .id <- row.names(model[[i]]$coord)

      # scores on the PC
      scores <- as.data.frame(model[[i]]$coord[,PC])
      names(scores) <- paste(".PC", PC, sep="")

      # square cosine : quality of the representation on the current space
      .cos2 <- model[[i]]$cos2[,PC]
      if (length(PC > 1)) {
        .cos2 <- rowSums(.cos2)
      }

      # contribution to the current PCs (for active variables and individuals only)
      if ("contrib" %in% names(model[[i]])) {
        .contrib <- model[[i]]$contrib[,PC]
        if (length(PC > 1)) {
          .contrib <- apply(.contrib, 1, function(x,v) {sum(x*v)}, explainedVar)
        }
      } else {
        .contrib <- 0
      }

      # TODO what is "cor" and what to do with it?

      # prepare result
      cD <- data.frame(.id, scores, .cos2, .contrib, .kind=i, stringsAsFactors=FALSE)

    } else {

      cD <- NULL

    }

    # store the result
    d <- c(d, list(cD))

  }

  # put everything in the same data.frame
  d <- do.call(rbind, d)


  # Add original data
  if (type == "observations") {
    # fetch data from the model object if necessary
    if (is.null(data)) {
      data <- model$call$X
    }
    # NB: we have to use join here and not cbind, in case there are supplementary observations
    data$.id <- row.names(data)
    d <- join(data, d, by=".id", type="full")
  }


  # store variance explained are attribute
  attr(d, "explained.variance") <- explainedVar

  return(d)
}

# TODO add a method for ade4

#' @method fortify pcaRes
#' @rdname fortify_pca
#' @export
fortify.pcaRes <- function(model, data=NULL, type=c("observations", "variables"), PC=c(1,2), ...) {
  #
  # Method for pcaMethods::pca
  #

  # Checks
  if (any(PC > ncol(model@scores))) {
    stop("At least one of the principal components does not exist")
  }
  if (length(PC) < 1) {
    stop("You must choose at least one principal component")
  }
  if (length(PC) > 2) {
    warning("Extracting information for more than two principal components. The plot might be difficult to read.")
  }
  type <- match.arg(type)


  # Compute variance explained
  # eigenvalues (scaled to be homogeneous with other packages)
  eig <- model@sDev^2
  # variance explained by each PC
  explainedVar <- eig / sum(eig)
  explainedVar <- explainedVar[PC]


  # Extract appropriate data
  if (type == "variables") {

    # Variables
    # variable identifier
    .id <- row.names(model@loadings)

    # loadings on the PCs
    loadings <- as.data.frame((t(t(model@loadings)*model@sDev)))
    # TODO select PCs here rather than later and save a bit of memory/cpu?
    # NB: loadings are scaled by std deviation, to be consistent with FactoMineR and ADE4, but this will result in a plot different from the one produced by stat::biplot, where loadings are not scaled
    names(loadings) <- paste(".", names(loadings), sep="")

    # squared cosine: quality of the representation on the current space
    .cos2 <- ( loadings / sqrt(rowSums(loadings^2)) )^2
    .cos2 <- .cos2[,PC]

    # contribution to the current PCs
    .contrib <- as.data.frame(t(t(loadings^2) / eig)) * 100
    .contrib <- .contrib[,PC]

    if (length(PC > 1)) {
      # the squared cos are additive
      .cos2 <- rowSums(.cos2)
      # contributions are scaled by the percentage of variance explained by each PC so that the contribution displayed is the contribution to the total variance projectable in the current space
      .contrib <- apply(.contrib, 1, function(x,v) {sum(x*v)}, explainedVar)
    }

    res <- data.frame(.id, loadings[,PC], .cos2, .contrib, .kind=type)

  }
  else if (type == "observations") {

    # Observations
    # observation identifier
    .id <- row.names(model@scores)
    if (is.null(.id)) {
      .id <- 1:nrow(model@scores)
    }

    # scores (i.e. coordinates) on the PCs
    scores <- as.data.frame(model@scores[,PC])
    names(scores) <- paste(".", names(scores), sep="")

    # square cosine
    .cos2 <- ( model@scores^2 / rowSums(model@scores^2) )
    .cos2 <- .cos2[,PC]

    # contributions
    .contrib <- as.data.frame( t( t(model@scores^2) * (1/nrow(model@scores)) / eig ) ) * 100
    .contrib <- .contrib[,PC]

    if (length(PC > 1)) {
      .cos2 <- rowSums(.cos2)
      .contrib <- apply(.contrib, 1, function(x,v) {sum(x*v)}, explainedVar)
    }

    res <- data.frame(.id, scores, .cos2, .contrib, .kind=type, stringsAsFactors=FALSE)

    # concatenate with original data
    if (is.null(data)) {
      data <- as.data.frame(model@completeObs)
    }
    # data can still be NULL because the @completeObs slot can be NULL in model
    if (!is.null(data)) {
      data$.id <- row.names(data)
      res <- join(data, res, by=".id", type="full")
      # TODO use the data in the object when it exists, for common columns, because it may have been imputed
    }
  }


  # store variance explained as attribute
  attr(res, "explained.variance") <- explainedVar

  return(res)
}

#' @method fortify rda
#' @rdname fortify_pca
#' @export
fortify.rda <- function(model, data=NULL, type=c("observations", "variables"), PC=c(1,2), ...) {
  #
  # Method for vegan::rda
  #
  # NB: those objects do not contain the original data, so data is NULL by default
  #

  # Checks
  if (any(PC > model$CA$rank)) {
    stop("At least one of the principal components does not exist")
  }
  if (length(PC) < 1) {
    stop("You must choose at least one principal component")
  }
  if (length(PC) > 2) {
    warning("Extracting information for more than two principal components. The plot might be difficult to read.")
  }
  type <- match.arg(type)


  # Compute variance explained
  # eigenvalues (scaled to be homogeneous with other packages)
  eig <- model$CA$eig
  # variance explained by each PC
  explainedVar <- eig / sum(eig)
  explainedVar <- explainedVar[PC]


  # Extract appropriate data
  if (type == "variables") {

    # Variables
    # variable identifier
    .id <- row.names(model$CA$v)

    # loadings on the PCs
    loadings <- as.data.frame((t(t(model$CA$v)*sqrt(model$CA$eig))))
    # TODO select PCs here rather than later and save a bit of memory/cpu?
    # NB: loadings are scaled by std deviation, to be consistent with FactoMineR and ADE4, but this will result in a plot different from the one produced by stat::biplot, where loadings are not scaled
    names(loadings) <- paste(".", names(loadings), sep="")

    # squared cosine: quality of the representation on the current space
    .cos2 <- ( loadings / sqrt(rowSums(loadings^2)) )^2
    .cos2 <- .cos2[,PC]

    # contribution to the current PCs
    .contrib <- as.data.frame(t(t(loadings^2) / eig)) * 100
    .contrib <- .contrib[,PC]

    if (length(PC > 1)) {
      # the squared cos are additive
      .cos2 <- rowSums(.cos2)
      # contributions are scaled by the percentage of variance explained by each PC so that the contribution displayed is the contribution to the total variance projectable in the current space
      .contrib <- apply(.contrib, 1, function(x,v) {sum(x*v)}, explainedVar)
    }

    res <- data.frame(.id, loadings[,PC], .cos2, .contrib, .kind=type)

  }
  else if (type == "observations") {

    # Observations
    # observation identifier
    .id <- row.names(model$CA$u)
    if (is.null(.id)) {
      .id <- 1:nrow(model$CA$u)
    }

    # scores (i.e. coordinates) on the PCs
    scores <- vegan::scores(model, choices=c(1:model$CA$rank), scaling=1)$sites
    # NB: requiring vegan::scores.rda should not be a problem since one needs vegan to produce the original object anyway
    # remove the scaling constant added by vegan, to make it compatible with prcomp
    scores <- scores * ( (nrow(scores)-1) * sum(model$CA$eig) ) ^ (1/4)
    colnames(scores) <- paste(".", colnames(scores), sep="")

    # square cosine
    .cos2 <- ( scores^2 / rowSums(scores^2) )
    .cos2 <- .cos2[,PC]

    # contributions
    .contrib <- as.data.frame( t( t(scores^2) * (1/nrow(scores)) / eig ) ) * 100
    .contrib <- .contrib[,PC]

    if (length(PC > 1)) {
      .cos2 <- rowSums(.cos2)
      .contrib <- apply(.contrib, 1, function(x,v) {sum(x*v)}, explainedVar)
    }

    res <- data.frame(.id, scores[,PC], .cos2, .contrib, .kind=type, stringsAsFactors=FALSE)
    if (!is.null(data)) {
      data$.id <- row.names(data)
      res <- join(data, res, by=".id", type="full")
    }
  }


  # store variance explained as attribute
  attr(res, "explained.variance") <- explainedVar

  return(res)
}
