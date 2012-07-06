\name{fortify.prcomp}
\alias{fortify.PCA}
\alias{fortify.pca}
\alias{fortify.pcaRes}
\alias{fortify.prcomp}
\title{Extract information from a Principal Component Analysis into a data.frame}
\usage{
  \method{fortify}{prcomp} (model, data = NULL,
    type = c("observations", "variables"), PC = c(1, 2),
    ...)

  \method{fortify}{PCA} (model, data = model$call$X,
    type = c("observations", "variables"), PC = c(1, 2),
    ...)

  \method{fortify}{pcaRes} (model,
    data = model@completeObs,
    type = c("observations", "variables"), PC = c(1, 2),
    ...)
}
\arguments{
  \item{model}{an object resulting from a PCA, using
  \code{\link[stats:prcomp]{stats::prcomp}},
  \code{\link[FactoMineR:PCA]{FactoMineR::PCA}} or
  \code{\link[pcaMethods:pca]{pcaMethods::pca}} (from
  bioconductor)}

  \item{data}{the original data used to compute the PCA, to
  be concatenated to the output when extracting
  observations; by default, the data will be extracted from
  the PCA object when possible (not for
  \code{\link[stats:prcomp]{prcomp}})}

  \item{type}{the type of data to extract : observations
  (i.e. rows, individuals) or variables (i.e. columns,
  descriptors); can be abbreviated}

  \item{PC}{the principal components to extract; more than
  two can be extracted but two is usual for plotting}

  \item{...}{pass-through argument}
}
\value{
  A data.frame containing the original data when pertinent
  and possible and the additional elements \item{.id}{the
  identifier of the observation (row name) or variable
  (column name)} \item{.PC#}{the scores (for observations)
  or loadings (for variables) on the extracted principal
  components} \item{.cos2}{the squared cosine summed over
  all extracted PCs} \item{.contrib}{the contribution to
  the selected PCs} \item{.kind}{the nature of the data
  extracted : observations, variables and possibly their
  status (active or supplementary)}
}
\description{
  Extract information from a Principal Component Analysis
  into a data.frame
}
\examples{
# PCA with stats::prcomp
pca <- prcomp(USArrests, scale = TRUE)

library("ggplot2")
head(fortify(pca))
fortify(pca, type = "variables")

# data is not containe in the `prcomp` object but can be provided
head(fortify(pca, data = USArrests, type = "observations"))

# select different principal components
fortify(pca, type = "var", PC=c(1,3))

\dontrun{
# PCA with FactoMineR::PCA
library("FactoMineR")

# add a missing value
d <- USArrests
d[1,2] <- NA
# use supplementary observations and variables
pca <- PCA(d, scale = TRUE, graph = FALSE, ind.sup = 2, quanti.sup = 4)

# the missing value is replaced by the column mean
# the supplementary observation is identified as such
head(fortify(pca))
head(fortify(pca, data = d))

# the supplementary variable is identified as such
fortify(pca, type = "variables")

# PCA with pcaMethods::pca, from bioconductor
library("pcaMethods")
pca <- pca(d, method = "svd", scale = "uv", completeObs = TRUE)

# the missing value is imputed by iterative PCA
head(fortify(pca))

fortify(pca, type="var")

}
}
\author{
  Jean-Olivier Irisson \email{irisson@normalesup.org}
}
\seealso{
  \code{\link[stats:prcomp]{prcomp}},
  \code{\link[FactoMineR:PCA]{PCA}}
}

