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


print.ggplot_list <- function(x) {
  #
  # Print a list of ggplots
  # Ask for the next plot when plotting interactively
  #
  
  # print the first plot
  print(x[[1]])
  
  # if there are more than one
  if (length(x) > 1) {
    # ask for the next plots if the device is interactive
    devAsk <- devAskNewPage()
    if (dev.interactive() || names(dev.cur()) == "null device") {
      devAskNewPage(TRUE)
    }

    # print the plots
    lapply(x[-1], print)

    # reset asking to the previous behaviour
    devAskNewPage(devAsk)
  }
  
  return(invisible(x))
}

plot.ggplot_list <- print.ggplot_list
 