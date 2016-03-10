# Determine the type of element to extract from an ordination object, allowing the vocabulary from different packages to co-exist
match_rowcol <- function(type, ...) {
  # define synonyms
  row_types <- c("rows", "lines", "observations", "individuals", "sites")
  col_types <- c("columns", "variables", "descriptors", "species")

  # allow abbreviation
  type <- match.arg(type, c(row_types, col_types), several.ok=FALSE)

  # reduce synonyms
  if ( type %in% row_types ) {
    type <- "row"
  } else {
    type <- "col"
  }
  # type[type %in% row_types] <- "row"
  # type[type %in% col_types] <- "col"
  # type <- unique(type)

  return(type)
}
