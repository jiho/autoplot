#' Automatic ggplot for a Linear Model
#'
#' @param object an object resulting from a linear model fit (of class \code{lm})
#'
#' @param type 
#'
#'
#' @return A ggplot2 object defining the plot
#'
#' @author Jean-Olivier Irisson \email{irisson@@normalesup.org}
#'
#' @seealso \link{fortify.lm}, \code{\link[stats:lm]{lm}}, \code{\link[stats:glm]{glm}}
#'
#' @examples
#' m <- lm(Rape ~ UrbanPop, data=USArrests)
#'
#' data <- fortify(m)
#'
#' }
#'

#' @method autoplot lm
#' @export
autoplot.lm <- function(object, type, mapping=aes(), ...) {

  # extrct data from the object
  data <- fortify(object, ...)
  
  # determine the type of plot to produce
  choices <- c("residuals-fitted", "normal qq", "scale-location", "cook's distance", "residuals-leverage", "cook's-leverage")
  
  if (is.numeric(type)) {
    type <- choices[type]
  } else {
    type <- match.arg(type, choices, several.ok=TRUE)
  }
  
  # produce plots
  plots <- list()
  if ("residual-fitted" %in% type) {
    p <- autoplot_lm_resid_fitted(data, mapping=mapping)
    plots <- c(plots, list(p))
  }

  if ("normal qq" %in% type) {
    p <- autoplot_lm_qq(data, mapping=mapping)
    plots <- c(plots, list(p))
  }

  if ("scale-location" %in% type) {
    p <- autoplot_lm_scale_loc(data, mapping=mapping)
    plots <- c(plots, list(p))
  }

  if ("cook's distance" %in% type) {
    p <- autoplot_lm_cooksd(data, mapping=mapping)
    plots <- c(plots, list(p))
  }

  if ("residuals-leverage" %in% type) {
    p <- autoplot_lm_resid_hat(data, mapping=mapping)
    plots <- c(plots, list(p))
  }

  if ("cook's-leverage" %in% type) {
    p <- autoplot_lm_hat_cooksd(data, mapping=mapping)
    plots <- c(plots, list(p))
  }
  
  class(plots) <- "ggplot_list"

  return(plots)
}


autoplot_lm_resid_fitted <- function(data, mapping, ...) {
  mapping <- c(mapping, aes(x=.fitted, y=.resid))
  class(mapping) <- "uneval"
  
  ggplot(data) + geom_point(mapping=mapping, ...)
}

autoplot_lm_qq <- function(data, mapping, ...) {
  map <- c(mapping, aes(sample=.stdresid))
  class(map) <- "uneval"
  
  ggplot(data) + geom_abline(intercept=0, slope=1, linetype="dotted") + stat_qq(mapping=map, ...)
}

autoplot_lm_scale_loc <- function(data, mapping, ...) {
  map <- c(mapping, aes(x=.fitted, y=sqrt(abs(.stdresid))))
  class(map) <- "uneval"
  
  ggplot(data) + geom_point(mapping=map, ...) + geom_smooth(mapping=map, se=FALSE, ...)
}

autoplot_lm_cooksd <- function(data, mapping, ...) {
  map <- c(mapping, aes(x=row.names(data), y=.cooksd))
  class(map) <- "uneval"
  
  ggplot(data) + geom_bar(mapping=map, stat="identity", width=0.2) + opts(axis.text.x=theme_text(angle=-45, hjust=0, vjust=1, colour="grey50", size=0.8*12))
}

autoplot_lm_resid_hat <- function(data, mapping, ...) {
  map <- c(mapping, aes(x=.hat, y=.stdresid))
  class(map) <- "uneval"
  
  ggplot(data, mapping=map) + geom_point() + geom_smooth(se=F)
}

autoplot_lm_hat_cooksd <- function(data, mapping, ...) {
  map <- c(mapping, aes(x=.hat, y=.cooksd))
  class(map) <- "uneval"
  
  ggplot(data, mapping=map) + geom_point() + geom_smooth(se=F) +
  geom_abline(intercept=0, slope=0, linetype="dotted") +
  geom_abline(intercept=0, slope=1, linetype="dotted") +
  geom_abline(intercept=0, slope=2, linetype="dotted") +
  geom_abline(intercept=0, slope=3, linetype="dotted")
}
