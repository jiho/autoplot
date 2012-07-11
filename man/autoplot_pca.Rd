\name{autoplot_pca}
\alias{autoplot.PCA}
\alias{autoplot.pcaRes}
\alias{autoplot.prcomp}
\alias{autoplot_pca}
\title{Automatic ggplot for a Principal Component Analysis}
\usage{
  autoplot_pca(object,
    type = c("observations", "variables"), mapping = aes(),
    data = NULL, PC = c(1, 2), ...)

  \method{autoplot}{prcomp} (object, ...)

  \method{autoplot}{PCA} (object, ...)

  \method{autoplot}{pcaRes} (object, ...)
}
\arguments{
  \item{object}{an object resulting from a PCA, using
  \code{\link[stats:prcomp]{stats::prcomp}} or
  \code{\link[FactoMineR:PCA]{FactoMineR::PCA}}}

  \item{type}{whether to extract observations (i.e. rows,
  individuals) or variables (i.e. columns, descriptors);
  can be abbreviated}

  \item{mapping}{a call to aes() specifying additional
  mappings between variables and plot aesthetics; by
  default, positions in x and y are mapped to the scores or
  loadings on the principal components and colour is mapped
  to the nature of the data (active or supplementary) when
  relevant. See \link{fortify_pca} for a list of the other
  mappable variables returned by the \code{fortify}
  methods}

  \item{data}{the original data used to compute the PCA, to
  be concatenated to the output when extracting
  observations. This allows to map original data columns to
  aesthetics of the plot, even if those columns were not
  used in the PCA. When \code{NULL} the data used in the
  PCA will be extracted from the PCA object, when possible
  (not for \code{\link[stats:prcomp]{prcomp}})}

  \item{PC}{the principal components to extract; two are
  necessary to produce a plot}

  \item{...}{passed to the various geoms; can be used to
  \emph{set} further aesthetics}
}
\value{
  A ggplot2 object defining the plot
}
\description{
  Automatic ggplot for a Principal Component Analysis
}
\examples{
# PCA with stats::prcomp
pca <- prcomp(USArrests, scale = TRUE)

library("ggplot2")
autoplot(pca)

# add further aesthetic mappings
names(fortify(pca, type = "obs"))
autoplot(pca, type = "obs", mapping=aes(alpha=.cos2))
autoplot(pca, type = "obs", mapping=aes(alpha=.cos2, size=.contrib))
# including from the original data
autoplot(pca, type = "obs", mapping=aes(alpha=.cos2, size=Murder), data=USArrests)

# aesthetics can also be set
autoplot(pca, type = "obs", mapping=aes(alpha=.cos2), colour="red", size=3, shape=15)

\dontrun{
# PCA with FactoMineR::PCA
library("FactoMineR")
# use supplementary observations and variables
pca <- PCA(USArrests, scale = TRUE, graph=FALSE, ind.sup = 2, quanti.sup = 4)

# colour is mapped by default
autoplot(pca)

# but the mapping can be overridden by mapping another variable
autoplot(pca, type = "obs", mapping = aes(colour=.contrib))
# or setting the corresponding aesthetic
autoplot(pca, type = "obs", colour="black")

# additional mappings can be specified
autoplot(pca, type = "obs", mapping = aes(colour=.contrib, alpha=.cos2, shape=.kind))
# in particular, with FactoMineR, the data is present by default and can be mapped
names(fortify(pca, type = "obs"))
autoplot(pca, "obs", aes(alpha=.cos2, size=Murder))

# PCA with pcaMethods::pca, from bioconductor
library("pcaMethods")
# equivalent to `prcomp`
pca <- pca(USArrests, method="svd", scale="uv", completeObs=TRUE, nPcs=4)
autoplot(pca)
autoplot(pca, type="obs", mapping=aes(alpha=.cos2))

# other computation method specific to pcaMethods
pca <- pca(USArrests, method="ppca", scale="uv", completeObs=TRUE, nPcs=4)
autoplot(pca)

}
}
\author{
  Jean-Olivier Irisson \email{irisson@normalesup.org}
}
\seealso{
  \link{fortify_pca}, \code{\link[stats:prcomp]{prcomp}} in
  package stats, \code{\link[FactoMineR:PCA]{PCA}} in
  package FactoMineR, \code{\link[pcaMethods:pca]{pca}} in
  package pcaMethods
}

