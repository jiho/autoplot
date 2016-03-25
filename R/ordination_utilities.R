# Remove centering and scaling done by function scale()
unscale <- function(x, center=NULL, scale=NULL) {
  if (is.null(center)) {
    center <- attr(x, "scaled:center")
  }
  if (is.null(scale)) {
    scale <- attr(x, "scaled:scale")
  }
  # TODO handle cases when only centering or only scaling is done
  if (is.null(center) | is.null(scale)) {
    stop("Cannot extract center and/or scale from the data and none is provided.")
  }
  unscaledx <- scale(x, center=(-center/scale), scale=1/scale)
  attr(unscaledx, "scaled:center") <- attr(unscaledx, "scaled:scale") <- NULL
  return(unscaledx)
}


# Determine the type of scores to extract from an ordination object
# This defines synonyms to allow the vocabulary from different packages to co-exist
match_type <- function(type, ...) {
  # define synonyms
  row_types <- c("rows", "lines", "observations", "objects", "individuals", "sites")
  col_types <- c("columns", "variables", "descriptors", "species")

  # allow abbreviation
  type <- match.arg(type, c(row_types, col_types, ...), several.ok=FALSE)
  # TODO consider allowing several.ok=TRUE to be able to match several elements of a vector

  # reduce synonyms
  type[type %in% row_types] <- "row"
  type[type %in% col_types] <- "col"

  return(type)
}


# Determine the type scaling for scores extracted from an ordination object
# This uses the same synonyms
match_scaling <- function(scaling) {
  # otherwise select a type of scaling
  if (is.numeric(scaling)) {
    # scaling can be a number, in a way compatible with most of the code out there (vegan in particular)
    if (! scaling %in% 0:3) {
      stop("Scaling should be 0, 1, 2, or 3")
    }
    scaling <- c("none", "row", "col", "both")[scaling+1]

  } else if (is.character(scaling)) {
    # or a character string, specifying in which "direction" to scale
    scaling <- match_type(scaling, c("none", "both", "biplot"))
    if (scaling == "biplot") { scaling <- "both" }

  } else {
    stop("Scaling should be a number between 0 and 3 or a character string")
  }
  
  return(scaling)
}
# TODO generalise to matches in a list cf choose plot

