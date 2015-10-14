#' @title Automatic ggplot for a Correspondence Analysis
#'
#' @param object an object resulting from a CA
#'
#' @param mapping a call to aes() specifying additional mappings between variables and plot aesthetics; by default, positions in x and y are mapped to the scores on the principal components and colour is mapped to the nature of the data (active or supplementary) when relevant. See \code{\link{fortify_ca}} for a list of the other mappable variables returned by the \code{fortify} methods
#'
#' @template ca_params
#'
#' @param ... passed to the various geoms; can be used to \emph{set} further aesthetics
#'
#' @return A ggplot2 object defining the plot
#'
#' @author Jean-Olivier Irisson \email{irisson@@normalesup.org}
#'
#' @seealso \code{\link{fortify_ca}} for the function which actually extracts the data from the CA object.
#'@name autoplot_ca
NULL

#' @method autoplot correspondence
#' @rdname autoplot_ca
#' @export
autoplot.correspondence <- function(ca, data=NULL, PC=c(1,2), mapping=aes(), ...) {

	# Extract data
	data = fortify.correspondence(ca, data=data, PC=PC)
	# NB: from now on, data is different from the data argument. Not very clean but useful
	
	# compute variance explained by each principal component
	eig = ca$cor^2
	explainedVar = eig / sum(eig)

	# Prepare legends and labels
	PCs = grep("PC", names(data), value=TRUE)
	axesLabels = paste(PCs, " (", format(explainedVar, digits=3), "%)", sep="")

	# compute full mapping from defaults + arguments
	map = c(aes_string(x=PCs[1], y=PCs[2], colour=".type"), mapping)
	class(map) = "uneval"

	# plot
	p = ggplot(data, mapping=map) +
		# points
		geom_point() +
		# point labels
		geom_text(aes(label=.id), size=3, vjust=-1) +
		# nice axes legends
		scale_x_continuous(axesLabels[1]) + scale_y_continuous(axesLabels[2])

	return(p)
}

#' @method autoplot CA
#' @rdname autoplot_ca
#' @export
autoplot.CA <- function(ca, data=NULL, PC=c(1,2), mapping=aes(), ...) {

	# Extract data
	data = fortify.CA(ca, data=data, PC=PC)
	# compute variance explained by each principal component
	explainedVar = ca$eig$`percentage of variance`[PC]

	# Prepare legends and labels
	PCs = grep("PC", names(data), value=TRUE)
	axesLabels = paste(PCs, " (", format(explainedVar, digits=3), "%)", sep="")

	# compute full mapping from defaults + arguments
	map = c(aes_string(x=PCs[1], y=PCs[2], colour=".type"), mapping)
	class(map) = "uneval"

	# plot
	p = ggplot(data, mapping=map) +
		# points
		geom_point() +
		# point labels
		geom_text(aes(label=.id), size=3, vjust=-1) +
		# nice axes legends
		scale_x_continuous(axesLabels[1]) + scale_y_continuous(axesLabels[2])

	return(p)
}
