#' Automatic ggplot for a Principal Component Analysis
#'
#' @param object an object resulting from a PCA, using \code{\link[stats:prcomp]{stats::prcomp}} or \code{\link[FactoMineR:PCA]{FactoMineR::PCA}}
#'
#' @param type whether to extract observations (i.e. rows, individuals) or variables (i.e. columns, descriptors); can be abbreviated
#'
#' @param mapping a call to aes() specifying additional mappings between variables and plot aesthetics; by default, positions in x and y are mapped to the scores or loadings on the principal components and colour is mapped to the nature of the data (active or supplementary) when relevant. See \link{fortify_pca} for a list of the other mappable variables returned by the \code{fortify} methods
#'
#' @param ... passed to \code{\link{fortify_pca}}, in particular to provide the original data if it cannot be retrieved from the object or to select wich principal components to plot
#'
#' @return A ggplot2 object defining the plot
#'
#' @author Jean-Olivier Irisson \email{irisson@@normalesup.org}
#'
#' @seealso \link{fortify_pca}, \code{\link[stats:prcomp]{prcomp}} in package stats, \code{\link[FactoMineR:PCA]{PCA}} in package FactoMineR, \code{\link[pcaMethods:pca]{pca}} in package pcaMethods
#'
#' @examples
#' # PCA with stats::prcomp
#' pca <- prcomp(USArrests, scale = TRUE)
#'
#' library("ggplot2")
#' autoplot(pca)
#'
#' # add further aesthetic mappings
#' names(fortify(pca, type = "obs"))
#' autoplot(pca, type = "obs", mapping=aes(alpha=.cos2))
#' autoplot(pca, type = "obs", mapping=aes(alpha=.cos2, size=.contrib))
#' # including from the original data
#' autoplot(pca, type = "obs", mapping=aes(alpha=.cos2, size=Murder), data=USArrests)
#'
#' \dontrun{
#' # PCA with FactoMineR::PCA
#' library("FactoMineR")
#' # use supplementary observations and variables
#' pca <- PCA(USArrests, scale = TRUE, graph=FALSE, ind.sup = 2, quanti.sup = 4)
#'
#' # colour is mapped by default
#' autoplot(pca)
#'
#' # but the mapping can be overridden
#' autoplot(pca, type = "obs", mapping = aes(colour=.contrib))
#' # and additional mappings can be specified
#' autoplot(pca, type = "obs", mapping = aes(colour=.contrib, alpha=.cos2, shape=.kind))
#'
#' # with FactoMineR, the data is present by default and can be mapped
#' names(fortify(pca, type = "obs"))
#' autoplot(pca, "obs", aes(alpha=.cos2, size=Murder))
#'
#' # PCA with pcaMethods::pca, from bioconductor
#' library("pcaMethods")
#' # equivalent to `prcomp`
#' pca <- pca(USArrests, method="svd", scale="uv", completeObs=TRUE, nPcs=4)
#' autoplot(pca)
#' autoplot(pca, type="obs", mapping=aes(alpha=.cos2))
#'
#' # other computation method specific to pcaMethods
#' pca <- pca(USArrests, method="ppca", scale="uv", completeObs=TRUE, nPcs=4)
#' autoplot(pca)
#'
#' }
#'
# TODO actually describe the plots, layers, mappings, as recommended on https://github.com/hadley/ggplot2/wiki/autoplot

#' @export
autoplot_pca <- function(object, type=c("observations", "variables"), mapping=aes(), ...) {

  # check arguments
  type <- choose_plots(type, choices=c("observations", "variables"))

  # prepare the appropriate plots
  p <- list()

  if ("observations" %in% type) {
    data <- fortify(object, type="observations", ...)
    p <- c(p, list(observations=autoplot_pca_obs(data=data, mapping=mapping)))
  }

  if ("variables" %in% type) {
    data <- fortify(object, type="variables", ...)
    p <- c(p, list(variables=autoplot_pca_vars(data=data, mapping=mapping)))
  }

  # give the output a special class with an appropriate print method
  class(p) <- c("ggplot_list", "ggplot", "list")

  return(p)
}

#' @method autoplot prcomp
#' @rdname autoplot_pca
#' @export
autoplot.prcomp <- function(object, type=c("observations", "variables"), mapping=aes(), ...) {
  autoplot_pca(object=object, type=type, mapping=mapping, ...)
}

#' @method autoplot PCA
#' @rdname autoplot_pca
#' @export
autoplot.PCA <- function(object, type=c("observations", "variables"), mapping=aes(), ...) {
  autoplot_pca(object=object, type=type, mapping=mapping, ...)
}

# TODO add a method for ade4

#' @method autoplot pcaRes
#' @rdname autoplot_pca
#' @export
autoplot.pcaRes <- function(object, type=c("observations", "variables"), mapping=aes(), ...) {
  autoplot_pca(object=object, type=type, mapping=mapping, ...)
}

autoplot_pca_axes_labels <- function(data) {
  # Prepare axes labels

  # get variance explained
  explVar <- attr(data, "explained.variance") * 100

  # get PC numbers
  PCs <- grep("PC", names(data), value=TRUE)
  PCs <- sub(".", "", PCs, fixed=TRUE)

  # concatenate with variance explained
  axesLabels <- paste(PCs, " (", format(explVar, digits=3), "%)", sep="")

  return(axesLabels)
}

autoplot_pca_vars <- function(data, mapping) {

  # Construct default aesthetic mappings
  # get PC numbers
  PCs <- grep("PC", names(data), value=TRUE)
  if (length(unique(data$.kind)) > 1) {
    # map colour to variable type (active or supplementary)
    mapping <- c(mapping, aes(colour=.kind))
    class(mapping) <- "uneval"
  }

  # Construct plot
  p <- ggplot(data, mapping=mapping)

  # set plot aspect
  p <- p +
    # square
    coord_fixed() +
    # plot a circle of radius 1
    geom_path(aes(x=cos(theta), y=sin(theta)), data=data.frame(theta=seq(0, 2*pi, length=100)), colour="white", size=0.5, linetype="solid", alpha=1) +
    # TODO adapt this to the the theme (currenty workd for theme_gray only)
    # TODO improve flexibility here: we are forced to set every unused aesthetic to avoid conflicts with the mappings inherited from the autoplot call
    # TODO R CMD check issues a note about theta being a global variable, this is probably because theta is never assigned outside of being a column of the data.frame and this caused the parser to chocke. setting "theta <- NULL" (or anything else) fixes the problem but is not very elegant
    xlim(-1.1,1.1) + ylim(-1.1,1.1)

  # plot data
  p <- p +
    # arrows describing the variables
    geom_segment(aes_string(x="0", y="0", xend=PCs[1], yend=PCs[2]), arrow=grid::arrow(angle=20, length=grid::unit(0.02, "npc"))) +
    # add variable names
    geom_text(aes_string(x=paste("1.04*", PCs[1], sep=""), y=paste("1.04*", PCs[2], sep=""), label=".id", hjust=paste("0.5-0.5*", PCs[1], sep=""), vjust=paste("0.5-0.5*", PCs[2], sep="")), size=3)
    # NB: the complex computation is to place the labels intelligently at the tip of the arrows

  # nice axes labels
  axesLabels <- autoplot_pca_axes_labels(data)
  p <- p + scale_x_continuous(axesLabels[1], breaks=seq(-1,1,0.5)) + scale_y_continuous(axesLabels[2], breaks=seq(-1,1,0.5))

  return(p)
}

autoplot_pca_obs <- function(data, mapping) {

  # Construct default aesthetic mappings
  # get PC numbers
  PCs <- grep("PC", names(data), value=TRUE)
  # map x/y position to PC
  mapping <- c(mapping, aes_string(x=PCs[1], y=PCs[2]))
  # map colour to variable type (active or supplementary) when necessary
  if (length(unique(data$.kind)) > 1) {
    mapping <- c(mapping, aes(colour=.kind))
  }
  class(mapping) <- "uneval"

  # Construct plot
  p <- ggplot(data, mapping=mapping)

  # plot data
  p <- p +
    # observations
    geom_point() +
    # labels
    geom_text(aes(label=.id), size=3, vjust=-1)

  # nice axes labels
  axesLabels <- autoplot_pca_axes_labels(data)
  p <- p + scale_x_continuous(axesLabels[1], breaks=0) + scale_y_continuous(axesLabels[2], breaks=0)

  return(p)
}
