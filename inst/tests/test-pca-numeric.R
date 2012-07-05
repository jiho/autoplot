context("Comparability of PCA functions")

# create some data
set.seed(123)
d <- data.frame(a=runif(5), b=runif(5), c=runif(5))

pcaS <- prcomp(d, scale=TRUE)
obsS <- fortify(pcaS, type="observations")
varS <- fortify(pcaS, type="variables")

library("FactoMineR")
pcaF <- PCA(d, scale.unit=TRUE, graph=F)
obsF <- fortify(pcaF, type="observations")
varF <- fortify(pcaF, type="variables")

test_that("eigenvalues are equal", {
  eigS <- pcaS$sdev^2
  eigF <- pcaF$eig$eigenvalue

  expect_equal(eigS, eigF)
})
