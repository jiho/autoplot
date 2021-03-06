#' @title Extract information from a Correspondence Analysis
#'
#' Extract information from a Correspondence Analysis such as the scores on the principal components and some fit statistics into a data.frame.
#'
#' @param model an object resulting from a CA
#'
#' @template ca_params
#'
#' @param ... pass-through argument
#'
#' @return
#' A data.frame containing
#'   \item{.id}{the identifier of the rows (row name) or columns (column name)}
#'   \item{.PC#}{the scores of rows or columns on the extracted principal components}
#'   \item{.cos2}{the squared cosine summed over all extracted PCs}
#'   \item{.contrib}{the contribution to the selected PCs}
#'   \item{.type}{the nature of the data extracted : rows, columns and possibly their status (active or supplementary)}
#'
#' @author Jean-Olivier Irisson \email{irisson@@normalesup.org}
#'
#' @seealso \code{\link{autoplot_ca}} to produce plots based on the output of \code{fortify}.
#'
#' @name fortify_ca
NULL


#' @method fortify correspondence
#' @rdname fortify_ca
#' @export
fortify.correspondence <- function(model, data=NULL, PC=c(1,2), ...) {
	# checks
	if (any(PC > ncol(model$rscore))) {
		stop("At least one of the principal components does not exist")
	}
	if (length(PC) < 1) {
		stop("You must choose at least one principal component")
	}
	if (length(PC) > 2) {
		warning("Extracting information for more than two principal components. The plot might be difficult to read.")
	}
	if (is.null(data)) {
		warning("fortify() needs to be provided with the original data (using the data argument)\n  to compute the contributions of rows and columns to the principal components.\n  Here they will be NA.")
	}

	# eigenvalues (scaled to be homogeneous with other packages)
	eig <- model$cor^2
	# variance explained by each PC
	explainedVar <- eig / sum(eig)

	d <- data.frame()

	for (dim in c("rscore", "cscore")) {
		.id <- row.names(model[[dim]])

		scores <- as.data.frame(t(t(model[[dim]]) * model$cor))
		names(scores) <- paste(".PC", 1:ncol(scores), sep="")

		.cos2 <- scores^2 / rowSums(scores^2)
		if (length(PC > 1)) {
			.cos2 <- rowSums(.cos2)
		}

		if (is.null(data)) {
			.contrib <- NA
		} else {
			if (dim == "rscore") {
				scaling <- rowSums(data) / sum(data)
			} else {
				scaling <- colSums(data) / sum(data)
			}
			.contrib <- as.data.frame(t(t(scores^2*scaling)/eig) * 100)
			if (length(PC > 1)) {
				.contrib <- apply(.contrib, 1, function(x,v) {sum(x*v)}, explainedVar[PC])
			}
		}
		
		.type <- ifelse(dim == "rscore", "row", "col")
		
		d <- rbind(d, data.frame(.id, scores[,PC], .cos2, .contrib, .type, stringsAsFactors=FALSE))
	}
	
	return(d)
}

#' @method fortify CA
#' @rdname fortify_ca
#' @export
fortify.CA <- function(model, data=NULL, PC=c(1,2), ...) {
	# checks
	if (any(PC > ncol(model$row$coord))) {
		stop("At least one of the principal components does not exist")
	}
	if (length(PC) < 1) {
		stop("You must choose at least one principal component")
	}
	if (length(PC) > 2) {
		warning("Extracting information for more than two principal components. The plot might be difficult to read.")
	}

	# variance explained by each PC (scaled to 1)
	explainedVar <- model$eig$eigenvalue / sum(model$eig$eigenvalue)

	d <- data.frame()

	for (i in c("row", "row.sup", "col", "col.sup")) {

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

			# contribution to the current PCs (for active variable only)
			if ("contrib" %in% names(model[[i]])) {
				.contrib <- model[[i]]$contrib[,PC]
				if (length(PC > 1)) {
					.contrib <- apply(.contrib, 1, function(x,v) {sum(x*v)}, explainedVar[PC])
				}
			} else {
				.contrib <- 0
			}

			cD <- data.frame(.id, scores, .cos2, .contrib, .type=i, stringsAsFactors=FALSE)

		} else {

			cD <- NULL

		}

		d <- rbind(d, cD)

	}

	return(d)
}
