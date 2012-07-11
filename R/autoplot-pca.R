#' Automatic ggplot for a Principal Component Analysis
#'
#' @param object an object resulting from a PCA, using \code{\link{prcomp}} in package \code{stats}, \code{\link{PCA}} in package \code{factoMineR} or \code{\link{pca}} in package \code{pcaMethods} from bioconductor
#'
#' @param type whether to extract observations (i.e. rows, individuals) or variables (i.e. columns, descriptors); can be abbreviated
#'
#' @param mapping a call to aes() specifying additional mappings between variables and plot aesthetics; by default, positions in x and y are mapped to the scores or loadings on the principal components and colour is mapped to the nature of the data (active or supplementary) when relevant. See \link{fortify_pca} for a list of the other mappable variables returned by the \code{fortify} methods
#'
#' @param data the original data used to compute the PCA, to be concatenated to the output when extracting observations. This allows to map original data columns to aesthetics of the plot, even if those columns were not used in the PCA. When \code{NULL} the data used in the PCA will be extracted from the PCA object, when possible (not for \code{\link{prcomp}})
#'
#' @param PC the principal components to extract; two are necessary to produce a plot
#'
#' @param ... passed to the various geoms; can be used to \emph{set} further aesthetics
#'
#' @return A ggplot2 object defining the plot
#'
#' @author Jean-Olivier Irisson \email{irisson@@normalesup.org}
#'
#' @seealso \link{fortify_pca}, \code{\link{prcomp}} in package \code{stats}, \code{\link{PCA}} in package \code{factoMineR} or \code{\link{pca}} in package \code{pcaMethods} from bioconductor
#'
#' @examples
#' # PCA with stats::prcomp
#' pca <- prcomp(USArrests, scale = TRUE)
#' autoplot(pca)
#'
#' # add further aesthetic mappings
#' names(fortify(pca, type = "obs"))
#' autoplot(pca, type = "obs", mapping=aes(alpha=.cos2))
#' autoplot(pca, type = "obs", mapping=aes(alpha=.cos2, size=.contrib))
#' # including from the original data
#' autoplot(pca, type = "obs", mapping=aes(alpha=.cos2, size=Murder), data=USArrests)
#'
#' # aesthetics can also be set
#' autoplot(pca, type = "obs", mapping=aes(alpha=.cos2), colour="red", size=3, shape=15)
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
#' # but the mapping can be overridden by mapping another variable
#' autoplot(pca, type = "obs", mapping = aes(colour=.contrib))
#' # or setting the corresponding aesthetic
#' autoplot(pca, type = "obs", colour="black")
#'
#' # additional mappings can be specified
#' autoplot(pca, type = "obs", mapping = aes(colour=.contrib, alpha=.cos2, shape=.kind))
#' # in particular, with FactoMineR, the data is present by default and can be mapped
#' names(fortify(pca, type = "obs"))
#' autoplot(pca, "obs", aes(alpha=.cos2, size=Murder))
#'
#' # PCA with pcaMethods::pca, from bioconductor
#' library("pcaMethods")
#' # settings equivalent `prcomp`
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
autoplot_pca <- function(object, type=c("observations", "variables"), mapping=aes(), data=NULL, PC=c(1, 2), ...) {

  # check arguments
  type <- choose_plots(type, choices=c("observations", "variables"))
  if (length(PC) != 2) {
    stop("You must choose exactly two principal components to plot")
  }

  # prepare the appropriate plots
  p <- list()

  if ("observations" %in% type) {
    fData <- fortify(model=object, data=data, type="observations", PC=PC)
    p <- c(p, list(observations=autoplot_pca_obs(data=fData, mapping=mapping, ...)))
  }

  if ("variables" %in% type) {
    fData <- fortify(model=object, data=data, type="variables", PC=PC)
    p <- c(p, list(variables=autoplot_pca_vars(data=fData, mapping=mapping, ...)))
  }

  if (length(p) == 1) {
    # when there is only one plot, just return it instead of a list
    # this allows to do
    #   autoplot() + geom_***
    # instead of having to do
    #   autoplot()[[1]] + geom_***
    # which is not very intuitive
    p <- p[[1]]
  } else {
    # give the output a special class with an appropriate print method
    class(p) <- c("ggplot_list", "ggplot", "list")
  }

  return(p)
}

#' @method autoplot prcomp
#' @rdname autoplot_pca
#' @export
autoplot.prcomp <- function(object, ...) {
  autoplot_pca(object=object, ...)
}

#' @method autoplot PCA
#' @rdname autoplot_pca
#' @export
autoplot.PCA <- function(object, ...) {
  autoplot_pca(object=object, ...)
}

# TODO add a method for ade4

#' @method autoplot pcaRes
#' @rdname autoplot_pca
#' @export
autoplot.pcaRes <- function(object, ...) {
  autoplot_pca(object=object, ...)
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

#' @importFrom grid arrow unit
autoplot_pca_vars <- function(data, mapping, ...) {

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
    geom_path(aes(x=cos(theta), y=sin(theta)), data=data.frame(theta=seq(0, 2*pi, length=100)), colour="white", size=0.5, linetype="solid", alpha=1)
    # TODO adapt this to the the theme (currenty workd for theme_gray only)
    # TODO improve flexibility here: we are forced to set every unused aesthetic to avoid conflicts with the mappings inherited from the autoplot call
    # TODO R CMD check issues a note about theta being a global variable, this is probably because theta is never assigned outside of being a column of the data.frame and this caused the parser to chocke. setting "theta <- NULL" (or anything else) fixes the problem but is not very elegant

  # plot data
  p <- p +
    # arrows describing the variables
    geom_segment(aes_string(x="0", y="0", xend=PCs[1], yend=PCs[2]), arrow=arrow(angle=20, length=unit(0.02, "npc")), ...) +
    # add variable names
    geom_text(aes_string(x=paste("1.04*", PCs[1], sep=""), y=paste("1.04*", PCs[2], sep=""), label=".id", hjust=paste("0.5-0.5*", PCs[1], sep=""), vjust=paste("0.5-0.5*", PCs[2], sep="")), ..., size=3)
    # NB: the complex computation is to place the labels intelligently at the tip of the arrows

  # nice axes labels
  axesLabels <- autoplot_pca_axes_labels(data)
  p <- p + scale_x_continuous(axesLabels[1], breaks=seq(-1,1,0.5), limits=c(-1.1,1.1), expand=c(0.1, 0)) + scale_y_continuous(axesLabels[2], breaks=seq(-1,1,0.5), limits=c(-1.1,1.1), expand=c(0.1, 0))

  return(p)
}

autoplot_pca_obs <- function(data, mapping, ...) {

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
    geom_point(...) +
    # labels
    geom_text(aes(label=.id), ..., size=3, vjust=-1)

  # nice axes labels
  axesLabels <- autoplot_pca_axes_labels(data)
  p <- p + scale_x_continuous(axesLabels[1], breaks=0, expand=c(0.1, 0)) + scale_y_continuous(axesLabels[2], breaks=0, expand=c(0.1, 0))

  return(p)
}
