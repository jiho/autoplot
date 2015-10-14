#' Automatic ggplot for a x, y, z list.
#'
#' Similar to \code{\link[graphics]{image}}
#'
#' @export
#' @param object list with components x, y and z
#' @param data not used by this method
#' @param mapping a call to aes() specifying additional mappings between variables and plot aesthetics; by default, positions in x and y are mapped to the x and y coordinates and fill is mapped to the values contained in the z component in the list
#' @param ... not used by this method
#' @examples
#' x <- seq(0, 1, by=0.1)
#' y <- seq(0, 2, by=0.1)
#' z <- outer(x, y, "+") 
#' d <- list(x=x, y=y, z=z)
#' image(d)
#' autoplot(d)
#' autoplot(d, mapping=aes(alpha=x))
autoplot.list <- function(object, data=NULL, mapping=aes(), ...) {
  d <- fortify(object)
  
  # add user mappings
  mapping <- c(mapping, aes_string(x="x", y="y", fill="z"))
  class(mapping) <- "uneval"
  
  p <- ggplot() + geom_tile(mapping=mapping, data=d)
  return(p)
}
