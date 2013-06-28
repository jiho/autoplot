
autoplot.correspondence <- function(ca, PC=c(1,2), data=NULL, mapping=aes(), ...) {

	# Extract data
	data = fortify.correspondence(ca, PC=PC, data=data)
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


autoplot.CA <- function(ca, PC=c(1,2), mapping=aes(), ...) {

	# Extract data
	data = fortify.CA(ca, PC=PC)
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


autoplot.MCA <- function(mca, PC=c(1,2), mapping=aes(), abbrev=FALSE, ...) {

	# Extract data
	data = fortify.MCA(mca, PC=PC)
	# compute variance explained by each principal component
	explainedVar = mca$eig$`percentage of variance`[PC]

	# Prepare legends and labels
	PCs = grep("PC", names(data), value=TRUE)
	axesLabels = paste(PCs, " (", format(explainedVar, digits=3), "%)", sep="")

	if (! abbrev) {
		# for qualitative variables, the levels became column header and this can be confusing if the names of the levels are not informative enough
		# so we give the option to have labels of the form
		#	name.level

		# original column names in the correct order
		qualiVar = mca$call$X[,c(mca$call$quali, mca$call$quali.sup)]
		qualiVarNames = names(qualiVar)
		# number of levels of each column
		qualiVarLevels = laply(qualiVar, nlevels)
		# create the vector of new names
		longIds = paste(rep(qualiVarNames, qualiVarLevels), data[data$.type %in% c("var", "quali.sup"), ".id"], sep=".")
		# store them as new Ids
		data[data$.type %in% c("var", "quali.sup"), ".id"] = longIds
	}

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
