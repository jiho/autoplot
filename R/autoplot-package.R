#' Automatic plotting with ggplot2
#'
#' ggplot2 is a versatile plotting package which allows to produce almost any kind of plot from data stored as a data.frame, by combining unit elements. However, it requires the user to design the plot entirely, from scratch. Many R functions for statistical analyses (linear models, factorial analyses, etc.) output objects of a given class and allow to easily plot classic diagnostics using \code{plot()}, by defining a specialized method for this generic function. This package aims at reproducing this functionality in ggplot2 while benefiting from its increased versatility. It provides two sets  of methods : (i) the \code{fortify()} methods extract data from the original object and format it as a data.frame; (ii) the \code{autoplot()} methods use these data.frames and leverage ggplot2 to produce diagnostic plots.
#'
#' @name autoplot-package
#' @docType package
#' @author Jean-Olivier Irisson \email{irisson@@normalesup.org}
#' @references
#' H. Wickham. Ggplot2: elegant graphics for data analysis. Springer-Verlag, New York, 2009.
#'
#' @keywords package
#'
#' @seealso See the original \code{\link[ggplot2:ggplot2-package]{ggplot2}} package.
#'
#' @importFrom ggplot2 autoplot fortify
NULL