#' Automatic ggplot for a Principal Component Analysis
#' 
#' Extract information from an object returned by a function performing Principal Component Analysis and produce a plot of observations, of variables, or a biplot.
#' 
#' @param object an object returned by a function performing Principal Component Analysis.
#'
#' @param mapping a call to aes() specifying additional mappings between variables and plot aesthetics. By default, positions \code{x} and \code{y} are mapped to the scores on the principal components.
#'
#' @param type the plot to produce: either plot "rows", "lines", "observations", "objects", "individuals", "sites" (which are all treated as synonyms), or plot "columns", "variables", "descriptors", "species" (which are, again, synonyms), or produce a "biplot". All can be abbreviated. By default, observations are plotted.
#'
#' @inheritParams pca_tidiers
#' @inheritParams scores
#'
#' @param n.max.labels maximum number of observation labels to plot. Let \code{n} be the number of data rows. For \code{0 < n < n.max.labels/2}, labels are plotted and shifted to avoid over-plotting. For \code{n.max.labels/2 < n < n.max.labels}, labels are plotted but not disambiguated. For \code{n.max.labels < n} labels are not plotted.
#'
#' @param ... passed to the various geoms; can be used to \emph{set} further aesthetics.
#'
#' @details
#' The object is passed to the appropriate \code{augment} method, defined in \code{\link{pca_tidiers}}, which extracts the scores, possibly the original data, and other relevant information from the PCA object. The resulting \code{data.frame} is plotted with:
#' \describe{
#'   \item{\code{\link[ggplot2]{geom_point}}:}{for observations points,}
#'   \item{\code{\link[ggplot2]{geom_text}} or \code{\link[ggrepel]{geom_text_repel}}:}{for observations labels,}
#'   \item{\code{\link[ggplot2]{geom_segment}}:}{for variables vectors.}
#' }
#' 
#' @return A ggplot2 object defining the plot.
#'
#' @template pca_seealso
#'
#' @examples
#' pca <- prcomp(USArrests, scale=TRUE)
#' autoplot(pca)
#' autoplot(pca, type="variables")
#' autoplot(pca, type="biplot") # defaults to scaling=3
#' autoplot(pca, type="biplot", scaling=1)
#' autoplot(pca, type="biplot", scaling=2, n.max.labels=0)
#'
#' # add further aesthetic mappings
#' names(augment(pca, data=USArrests))
#' autoplot(pca, data=USArrests, mapping=aes(alpha=.cos2))
#' autoplot(pca, data=USArrests, mapping=aes(alpha=.cos2, size=.contrib))
#' # including from the original data
#' autoplot(pca, data=USArrests, mapping=aes(alpha=.cos2, color=Murder))
#'
#' # aesthetics can also be set
#' autoplot(pca, mapping=aes(alpha=.cos2), color="darkblue", size=3)
#'
#' if (require("FactoMineR")) {
# TODO use supp variables
# #'   # use supplementary observations and variables
# #'   pca <- PCA(USArrests, scale = TRUE, graph=FALSE, ind.sup = 2, quanti.sup = 4)
#'   pca <- PCA(USArrests, graph=FALSE)
#'   autoplot(pca)
#'   autoplot(pca, type="variables")
#'   # with FactoMineR, the data is present by default and can be mapped
#'   names(augment(pca))
#' autoplot(pca, mapping=aes(alpha=.cos2, color=Murder))
# #'   # colour is mapped by default
# #'   
# #'   # but the mapping can be overridden by mapping another variable
# #'   autoplot(pca, type = "obs", mapping = aes(colour=.contrib))
# #'   # or setting the corresponding aesthetic
# #'   autoplot(pca, type = "obs", colour="black")
# #'   
# #'   # additional mappings can be specified
# #'   autoplot(pca, type = "obs", mapping = aes(colour=.contrib, alpha=.cos2, shape=.type))
#' }
#'
#' if (require("vegan")) {
#'   pca <- vegan::rda(USArrests, scale=TRUE)
#'   plot(pca)
#'   autoplot(pca, type="biplot", scaling=2)
#'   autoplot(pca)
#' }
#'
#' if (require("ade4")) {
#'   pca <- ade4::dudi.pca(USArrests, nf=4, scannf=FALSE)
#'   biplot(pca)
#'   autoplot(pca, type="biplot", scaling=1)
#' }
#'
#' if (require("pcaMethods")) {
#'   pca <- pca(USArrests, scale="uv", nPcs=4)
#'   biplot(pca)
#'   autoplot(pca, type="biplot")
#' }
#'
#'@name autoplot_pca
NULL

autoplot_pca <- function(object, mapping=aes(), data=NULL, dimensions=c(1,2), type="rows", scaling=type, n.max.labels=100, ...) {

  # TODO use which instead of type, see plot.lm
  # TODO actually describe the plots, layers, mappings, as recommended on https://github.com/hadley/ggplot2/wiki/autoplot

  # check arguments
  type <- match_type(type, c("biplot"))
  scaling <- match_scaling(scaling)

  dimensions <- na.omit(dimensions)
  if (length(dimensions) != 2) {
    stop("Argument 'dimensions' should be a vector of length two: the two principal components to plot.")
  }

  # prepare default plot
  p <- ggplot() + coord_fixed() + scale_x_continuous(breaks=0, expand=c(0.2, 0))  + scale_y_continuous(breaks=0, expand=c(0.2, 0)) 
  # and axes labels
  eig <- tidy(object)[dimensions,]
  axes_names <- ordination_axes_titles(eig)
  p <- p + labs(x=axes_names[1], y=axes_names[2])
  
  
  # prepare the appropriate plots
  if (type == "row") {
    # get data
    d <- augment(x=object, data=data, type="row", dimensions=dimensions, scaling=1)

    # build plot
    mapping <- ordination_mapping(data=d, mapping=mapping)
    p <- p + geom_ordination_points(data=d, mapping=mapping, n.max.labels=n.max.labels, ...)
  }

  if (type == "col") {
    d <- augment(x=object, data=data, type="col", dimensions=dimensions, scaling=0)

    # add the unit circle (and make it look similar to grid lines)
    # NB: getting the theme element works at the time the plot is created, so it will work for constructs such as
    #     theme_set(theme_bw())
    #     autoplot(...)
    # but not for
    #     autoplot(...) + theme_bw()
    # because at the time the autoplot function is called, the theme is not yet set to theme_bw()
    # TODO add circle for significance threshold
    grid <- calc_element("panel.grid.major", theme_get())
    p <- p + geom_path(aes_string(x="cos(theta)", y="sin(theta)"), data=data.frame(theta=seq(0, 2*pi, length=100)), colour=calc_element("panel.grid.major", theme_get())$colour, size=grid$size, linetype=grid$linetype, lineend=grid$lineend, alpha=1)
    # p <- p + coord_fixed(xlim=c(-1, 1), ylim=c(-1, 1))
    
    # build plot
    mapping <- ordination_mapping(data=d, mapping=mapping)
    p <- p + geom_ordination_vectors(data=d, mapping=mapping, ...)
    
  }

  if (type == "biplot") {
    dr <- augment(x=object, data=data, type="row", dimensions=dimensions, scaling=scaling)
    dc <- augment(x=object, data=data, type="col", dimensions=dimensions, scaling=scaling)
    # build plot
    mapping <- ordination_mapping(data=dr, mapping=mapping)
    p <- p +
      geom_ordination_points(data=dr, mapping=mapping, n.max.labels=n.max.labels, ...) +
      geom_ordination_vectors(data=dc, mapping=mapping, ...)
  }

  # if (length(p) == 1) {
  #   # when there is only one plot, just return it instead of a list
  #   # this allows to do
  #   #   autoplot() + geom_***
  #   # instead of having to do
  #   #   autoplot()[[1]] + geom_***
  #   # which is not very intuitive
  #   p <- p[[1]]
  # } else {
  #   # give the output a special class with an appropriate print method
  #   class(p) <- c("ggplot_list", "ggplot", "list")
  # }

  return(p)
}

ordination_mapping <- function(data, mapping, ...) {
  # get PC numbers
  PCs <- grep(".PC", names(data), value=TRUE)
  # map x/y position to PC and add other mappings
  mapping <- as.aes(mapping, aes_string(x=PCs[1], y=PCs[2]))
  # TODO deal with supplementary variables
  # map colour to variable type (active or supplementary) when necessary
  # if (length(unique(data$.type)) > 1) {
  #   mapping <- c(mapping, aes(colour=.type))
  # }
}

ordination_axes_titles <- function(eig) {
  paste0(eig$term, "(", round(eig$prop.variance * 100), "%)")
}

# Plot of points and labels in ordination space
geom_ordination_points <- function(data, mapping, n.max.labels, ...) {
  g <- geom_point(mapping=mapping, data=data, ...)
  
  # add data labels
  # few: try to avoid overlap with geom_text_repel
  # many: just plot text
  # too many: do not add labels altogether
  # TODO capture the size aes and rescale it to for text match a "natural" size for points
  n <- nrow(data)
  if (n <= (n.max.labels / 2)) {
    g <- c(
      g,
      ggrepel::geom_text_repel(as.aes(mapping, aes_string(label=".rownames")), data=data, ..., size=3)
    )
    # TODO capture color and map it to the segment color (using alpha)
  } else if (n <= n.max.labels) {
    g <- c(
      g,
      geom_text(as.aes(mapping, aes_string(label=".rownames")), data=data, ..., size=3, nudge_x=0.03, hjust=0)
    )
  }
  
  return(g)
}

# Plot of variables vectors in ordination space
geom_ordination_vectors <- function(data,  mapping, ...) {
  g <- list(
    # arrows describing the variables
    geom_segment(as.aes(aes_string(x="0", y="0", xend=mapping$x, yend=mapping$y), mapping), data=data, ..., lineend="round"),
    # add variable names
    geom_text(as.aes(aes_string(x=mapping$x, y=mapping$y, label=".rownames", hjust=paste0("ifelse(", mapping$x ,">0, -0.05, 1.05)"), mapping)), data=data, ..., size=3, vjust=0.5)
    # or use a "complex"" computation to place the labels intelligently at the tip of the arrows
    # geom_text(as.aes(aes_string(x=paste("1.04*", mapping$x, sep=""), y=paste("1.04*", mapping$y, sep=""), label=".label", hjust=paste("0.5-1*", mapping$x, sep=""), vjust=paste("0.5-1*", mapping$y, sep="")), mapping), data=data, ..., size=3)
  )
  return(g)
}

# Define the methods

#' @method autoplot prcomp
#' @rdname autoplot_pca
#' @export
autoplot.prcomp <- autoplot_pca

#' @method autoplot PCA
#' @rdname autoplot_pca
#' @export
autoplot.PCA <- autoplot_pca

#' @method autoplot rda
#' @rdname autoplot_pca
#' @export
autoplot.rda <- autoplot_pca

#' @method autoplot pca
#' @rdname autoplot_pca
#' @export
autoplot.pca <- autoplot_pca

#' @method autoplot pcaRes
#' @rdname autoplot_pca
#' @export
autoplot.pcaRes <- autoplot_pca
