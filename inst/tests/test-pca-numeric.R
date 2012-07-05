context("Comparability of PCA functions")

# create some data
set.seed(123)
d <- data.frame(a=runif(5), b=runif(5), c=runif(5))

# PCA with stats::prcomp
pcaS <- prcomp(d, scale=TRUE)
obsS <- fortify(pcaS, type="observations")
varS <- fortify(pcaS, type="variables")

# PCA with FactoMineR::PCA
library("FactoMineR")
pcaF <- PCA(d, scale.unit=TRUE, graph=F)
obsF <- fortify(pcaF, type="observations")
varF <- fortify(pcaF, type="variables")

# TODO add ade4

test_that("eigenvalues are equal", {
  eigS <- pcaS$sdev^2
  eigF <- pcaF$eig$eigenvalue

  expect_equal(eigS, eigF)
})

# TODO add some more tests, comparing the values returned by the various fortify methods. Ideally, they should all be equivalent.
