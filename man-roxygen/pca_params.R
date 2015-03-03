#'
#' @param data the original dataset, to be concatenated to the output when extracting observations. Whhen \code{NULL}, the default, data will be extracted from the PCA object when it contains it (i.e. not for \code{\link{prcomp}})
#'
#' @param type the elements to consider: either "observations", "individuals", "sites", "lines" (which are all synonyms) or "variables", "descriptors", "species", "columns" (which are, again, synonyms)
#'
#' @param PC the principal components to extract; two is usual for plotting
#'
#' @param scaling how to scale scores. Can be "none", "obs", "var", "both" (or synonyms of "obs" and "var", as for the argument \code{type}). Can also be a integer from 0 to 3, which designs a scaling as in the order above (0="none", 1="obs", etc.). Scaling observations makes the distances in the plot of observations approximations of the euclidan distance in the original space. Scaling variables relates the angle between variables to the correlations between them (right angle = no correlation). Scaling of observation or variables is by the square root of the eigenvalues. Scaling both is by the fourth root of the eigenvalues.
#' TODO add a note about auto scaling
#'
