\name{print.ggplot_list}
\alias{plot.ggplot_list}
\alias{print.ggplot_list}
\title{Draw a several plots successively on current graphics device}
\usage{
  \method{print}{ggplot_list} (x, ...)

  \method{plot}{ggplot_list} (x, ...)
}
\arguments{
  \item{x}{list of ggplots to display}

  \item{...}{other arguments not used by this method}
}
\description{
  When plotting several plots on an interactive device, the
  user is given the chance to examine each successive plot
  and has to press enter to see the next one. When plotting
  to a file, all plots are produced in sequence, without
  user interaction. This results in several pages in a PDF
  file. For image files (JPEG, PNG, etc.) the file name has
  to contain a numeric placeholder to produce several files
  instead of overwritting them (e.g. "Rplot%03d.png").
}

