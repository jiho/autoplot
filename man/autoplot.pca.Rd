\name{autoplot.prcomp}
\alias{autoplot.PCA}
\alias{autoplot.pca}
\alias{autoplot.pcaRes}
\alias{autoplot.prcomp}
\title{Automatic ggplot for a Principal Component Analysis}
\usage{
  \method{autoplot}{prcomp} (object,
    type = c("observations", "variables"), mapping = aes(),
    ...)

  \method{autoplot}{PCA} (object,
    type = c("observations", "variables"), mapping = aes(),
    ...)

  \method{autoplot}{pcaRes} (object,
    type = c("observations", "variables"), mapping = aes(),
    ...)
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
  relevant. See \link{fortify.pca} for a list of the other
  mappable variables returned by the \code{fortify}
  methods}

  \item{...}{passed to \code{\link{fortify.pca}}, in
  particular to provide the original data if it cannot be
  retrieved from the object or to select wich principal
  components to plot}
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

\dontrun{
# PCA with FactoMineR::PCA
library("FactoMineR")
# use supplementary observations and variables
pca <- PCA(USArrests, scale = TRUE, graph=FALSE, ind.sup = 2, quanti.sup = 4)

# colour is mapped by default
autoplot(pca)

# but the mapping can be overridden
autoplot(pca, type = "obs", mapping = aes(colour=.contrib))
# and additional mappings can be specified
autoplot(pca, type = "obs", mapping = aes(colour=.contrib, alpha=.cos2, shape=.kind))

# with FactoMineR, the data is present by default and can be mapped
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
  \link{fortify.pca}, \code{\link[stats:prcomp]{prcomp}},
  \code{\link[FactoMineR:PCA]{PCA}}
}

