choose_plots <- function(type, choices) {
  #
  # Expands the argument `type` into the values in choices
  # type can be numeric, abbreviations of choices, or full elements of choices
  # 

  # numeric = index in the choices vector
  if (is.numeric(type)) {
    # check the indexes are availables
    validTypes <- type %in% 1:length(choices)
    if (any(!validTypes)) {
      stop("Cannot find a plot type matching : ", paste(type[!validTypes], collapse=", "))
    }
    # extract the corresponding character representation
    type <- choices[type]
  }
  
  # character
  else if (is.character(type)){
    # allow abbreviations
    type <- match.arg(type, choices, several.ok=TRUE)
  } 

  else {
    stop("Cannot understand a type argument of class : ", class(type))
  }
  
  return(type)
}

#' Draw a several plots successively on current graphics device
#'
#' When plotting several plots on an interactive device, the user is given the chance to examine each successive plot and has to press enter to see the next one.
#' When plotting to a file, all plots are produced in sequence, without user interaction. This results in several pages in a PDF file. For image files (JPEG, PNG, etc.) the file name has to contain a numeric placeholder to produce several files instead of overwritting them (e.g. "Rplot%03d.png").
#'
#' @param x list of ggplots to display
#' @param ... other arguments not used by this method
#'
#' @method print ggplot_list
#' @export
print.ggplot_list <- function(x, ...) {
  #
  # Print a list of ggplots
  # Ask for the next plot when plotting interactively
  #
  
  # print the first plot
  print(x[[1]])
  
  # if there are more than one
  if (length(x) > 1) {

    # ask for the next plots if the device is interactive
    if (dev.interactive() || names(dev.cur()) == "null device") {
      # previous behaviour
      devAsk <- devAskNewPage()
      # turn asking on
      devAskNewPage(TRUE)
      # reset asking to the previous behaviour
      on.exit(devAskNewPage(devAsk))
    }

    # print the plots
    lapply(x[-1], print)

  }
  
  return(invisible(x))
}

#' @method plot ggplot_list
#' @rdname print.ggplot_list
#' @export
plot.ggplot_list <- print.ggplot_list
 