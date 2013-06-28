#
# Extract information from factorial analyses (PCA, CA, etc.)
#
# Different functions can perform these analyses in R. Several
# are supported here, and their results are homogenised.
#
# (c) Copyright 2011 Jean-Olivier Irisson. GNU General Public License
#
#------------------------------------------------------------


fortify.prcomp <- function(pca, PC=c(1,2))
#
#	Extract information from a PCA made using stats::prcomp
#
#	Arguments
#	PC		principal components to extract
#
#	Value
#	A list with two data.frames
#	variables      for the projection of variables
#	observations   for the projection of the observations
#
{
	# checks
	if (any(PC > ncol(pca$x))) {
		stop("At least one of the principal components does not exist")
	}
	if (length(PC) < 1) {
		stop("You must choose at least one principal component")
	}
	if (length(PC) > 2) {
		warning("Extracting information for more than two principal components. The plot might be difficult to read.")
	}

	# eigenvalues (scaled to be homogeneous with other packages)
	eig = pca$sdev^2
	# variance explained by each PC
	explainedVar = eig / sum(eig)

	# Variables
	# variable identifier
	.id = row.names(pca$rotation)

	# loadings on the PCs
	loadings = as.data.frame((t(t(pca$rotation)*pca$sdev)))
	# NB: loadings are scaled by std deviation, to be consistent with FactoMineR and ADE4, but this is different from the plot made by stat::biplot, where loadings are not scaled
	names(loadings) = paste(".", names(loadings), sep="")

	# squared cosine: quality of the representation on the current space
	.cos2 = ( loadings / sqrt(rowSums(loadings^2)) )^2
	.cos2 = .cos2[,PC]

	# contribution to the current PCs
	.contrib = as.data.frame(t(t(loadings^2) / eig)) * 100
	.contrib = .contrib[,PC]

	if (length(PC > 1)) {
		# the squared cos are additive
		.cos2 = rowSums(.cos2)
		# contributions are scaled by the percentage of variance explained by each PC so that the contribution displayed is the contribution to the total variance projectable in the current space
		.contrib = apply(.contrib, 1, function(x,v) {sum(x*v)}, explainedVar[PC])
	}

	variables = data.frame(.id, loadings[,PC], .cos2, .contrib, .type="var", stringsAsFactors=FALSE)

	# Observations
	# observation identifier
	.id = row.names(pca$x)
	# TODO: understand the difference between those x and the ind$coord of FactoMineR

	# scores on the PCs
	scores = as.data.frame(pca$x[,PC])
	names(scores) = paste(".", names(scores), sep="")

	# square cosine on the current space
	.cos2 = ( pca$x^2 / rowSums(pca$x^2) )
	.cos2 = .cos2[,PC]

	# contributions on the current space
	.contrib = as.data.frame( t( t(pca$x^2) * (1/nrow(pca$x)) / eig ) ) * 100
	.contrib = .contrib[,PC]

	if (length(PC > 1)) {
		.cos2 = rowSums(.cos2)
		.contrib = apply(.contrib, 1, function(x,v) {sum(x*v)}, explainedVar[PC])
	}

	observations = data.frame(.id, scores, .cos2, .contrib, .type="obs", stringsAsFactors=FALSE)
	# NB: objects of class prcomp do not contain the original data, so it cannot be included here, as is usually done in fortify.*

	return(list(variables=variables, observations=observations))
}


fortify.PCA <- function(pca, PC=c(1,2))
#
#	Extract information from a PCA made using FactoMineR::PCA
#
#	Arguments
#	PC		principal components to extract
#
#	Value
#	A list with two data.frames
#	variables      for the projection of variables
#	observations   for the projection of the observations + original data
#
{
  
	# checks
	if (any(PC > ncol(pca$ind$coord))) {
		stop("At least one of the principal components does not exist")
	}
	if (length(PC) < 1) {
		stop("You must choose at least one principal component")
	}
	if (length(PC) > 2) {
		warning("Extracting information for more than two principal components. The plot might be difficult to read.")
	}

	# variance explained by each PC (scaled to 1)
	explainedVar = pca$eig$eigenvalue / sum(pca$eig$eigenvalue)

	d <- list()

	for (i in c("var", "quanti.sup", "ind", "ind.sup", "quali.sup")) {

		if (!is.null(pca[[i]])) {

			# variable identifier
			.id = row.names(pca[[i]]$coord)

			# scores on the PC
			scores = as.data.frame(pca[[i]]$coord[,PC])
			names(scores) = paste(".PC", PC, sep="")

			# square cosine : quality of the representation on the current space
			.cos2=pca[[i]]$cos2[,PC]
			if (length(PC > 1)) {
				.cos2 = rowSums(.cos2)
			}

			# contribution to the current PCs (for active variables and individuals only)
			if ("contrib" %in% names(pca[[i]])) {
				.contrib=pca[[i]]$contrib[,PC]
				if (length(PC > 1)) {
					.contrib = apply(.contrib, 1, function(x,v) {sum(x*v)}, explainedVar[PC])
				}
			} else {
				.contrib = 0
			}

			# TODO what is "cor" and what to do with it?

			cD = data.frame(.id, scores, .cos2, .contrib, .type=i, stringsAsFactors=FALSE)

		} else {

			cD = NULL

		}

		d <- c(d, list(cD))

	}

	# variables + supplementary quantitative variables
	var = do.call(rbind, d[1:2])

	# individuals + supplementary individuals + supplementary qualitative variables
	obs = do.call(rbind, d[3:5])
	# append the original data
	X = pca$call$X
	X$.id = row.names(X)
	require("plyr")
	obs = join(obs, X, by=".id", type="full")

	return(list(variables=var, observations=obs))
}


fortify.correspondence <- function(ca, PC=c(1,2), data=NULL)
#
#	Extract information from a CA made using MASS::corresp
#
#	Arguments
#	PC		principal components to extract
#	data	the original data (which is necessary to compute the contributions)
#
#	Value
#	A data.frame with the projection of row and columns
#
{
	# checks
	if (any(PC > ncol(ca$rscore))) {
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
	eig = ca$cor^2
	# variance explained by each PC
	explainedVar = eig / sum(eig)

	d = data.frame()

	for (dim in c("rscore", "cscore")) {
		.id = row.names(ca[[dim]])

		scores = as.data.frame(t(t(ca[[dim]]) * ca$cor))
		names(scores) = paste(".PC", 1:ncol(scores), sep="")

		.cos2 = scores^2 / rowSums(scores^2)
		if (length(PC > 1)) {
			.cos2 = rowSums(.cos2)
		}

		if (is.null(data)) {
			.contrib = NA
		} else {
			if (dim == "rscore") {
				scaling = rowSums(data) / sum(data)
			} else {
				scaling = colSums(data) / sum(data)
			}
			.contrib = as.data.frame(t(t(scores^2*scaling)/eig) * 100)
			if (length(PC > 1)) {
				.contrib = apply(.contrib, 1, function(x,v) {sum(x*v)}, explainedVar[PC])
			}
		}
		
		.type = ifelse(dim == "rscore", "row", "col")
		
		d = rbind(d, data.frame(.id, scores[,PC], .cos2, .contrib, .type, stringsAsFactors=FALSE))
	}
	
	return(d)
}


fortify.CA <- function(ca, PC=c(1,2))
#
#	Extract information from a CA made using FactoMineR::CA
#
#	Arguments
#	PC		principal components to extract
#
#	Value
#	A data.frame with the projection of row and columns
#
{
	# checks
	if (any(PC > ncol(ca$row$coord))) {
		stop("At least one of the principal components does not exist")
	}
	if (length(PC) < 1) {
		stop("You must choose at least one principal component")
	}
	if (length(PC) > 2) {
		warning("Extracting information for more than two principal components. The plot might be difficult to read.")
	}

	# variance explained by each PC (scaled to 1)
	explainedVar = ca$eig$eigenvalue / sum(ca$eig$eigenvalue)

	d <- data.frame()

	for (i in c("row", "row.sup", "col", "col.sup")) {

		if (!is.null(ca[[i]])) {

			# variable identifier
			.id = row.names(ca[[i]]$coord)

			# scores on the PC
			scores = as.data.frame(ca[[i]]$coord[,PC])
			names(scores) = paste(".PC", PC, sep="")

			# square cosine : quality of the representation on the current space
			.cos2=ca[[i]]$cos2[,PC]
			if (length(PC > 1)) {
				.cos2 = rowSums(.cos2)
			}

			# contribution to the current PCs (for active variable only)
			if ("contrib" %in% names(ca[[i]])) {
				.contrib=ca[[i]]$contrib[,PC]
				if (length(PC > 1)) {
					.contrib = apply(.contrib, 1, function(x,v) {sum(x*v)}, explainedVar[PC])
				}
			} else {
				.contrib = 0
			}

			cD = data.frame(.id, scores, .cos2, .contrib, .type=i, stringsAsFactors=FALSE)

		} else {

			cD = NULL

		}

		d <- rbind(d, cD)

	}

	return(d)
}


fortify.MCA <- function(mca, PC=c(1,2))
#
#	Extract information from a MCA made using FactoMineR::MCA
#
#	Arguments
#	PC		principal components to extract
#
#	Value
#
{
	# checks
	if (any(PC > ncol(mca$ind$coord))) {
		stop("At least one of the principal components does not exist")
	}
	if (length(PC) < 1) {
		stop("You must choose at least one principal component")
	}
	if (length(PC) > 2) {
		warning("Extracting information for more than two principal components. The plot might be difficult to read.")
	}

	# variance explained by each PC (scaled to 1)
	explainedVar = mca$eig$eigenvalue / sum(mca$eig$eigenvalue)

	d <- data.frame()

	for (i in c("var", "quanti.sup", "ind", "ind.sup", "quali.sup")) {

		if (!is.null(mca[[i]])) {

			# variable identifier
			.id = row.names(mca[[i]]$coord)

			# scores on the PC
			scores = as.data.frame(mca[[i]]$coord[,PC])
			names(scores) = paste(".PC", PC, sep="")

			# square cosine : quality of the representation on the current space
			.cos2=mca[[i]]$cos2[,PC]
			if (length(PC > 1)) {
				.cos2 = rowSums(.cos2)
			}

			# contribution to the current PCs (for active variables and individuals only)
			if ("contrib" %in% names(mca[[i]])) {
				.contrib=mca[[i]]$contrib[,PC]
				if (length(PC > 1)) {
					.contrib = apply(.contrib, 1, function(x,v) {sum(x*v)}, explainedVar[PC])
				}
			} else {
				.contrib = 0
			}

			# TODO what to do with eta2 ?
			# TODO what to do with v.test ?

			cD = data.frame(.id, scores, .cos2, .contrib, .type=i, stringsAsFactors=FALSE)

		} else {

			cD = NULL

		}

		d <- rbind(d, cD)

	}

	# append the original data
	X = mca$call$X
	X$.id = row.names(X)
	require("plyr")
	d = join(d, X, by=".id", type="full")

	return(d)
}

