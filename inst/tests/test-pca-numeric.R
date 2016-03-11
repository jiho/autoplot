context("Comparability of PCA functions")

# create some data
set.seed(123)
d <- data.frame(a=runif(5), b=runif(5), c=runif(5))

# PCA with stats::prcomp
pcaS <- prcomp(d, scale=TRUE)
obsS <- fortify(pcaS, type="observations")
varS <- fortify(pcaS, type="variables")

# PCA with FactoMineR::PCA
suppressPackageStartupMessages(library("FactoMineR"))
pcaF <- PCA(d, scale.unit=TRUE, graph=F)
obsF <- fortify(pcaF, type="observations")
varF <- fortify(pcaF, type="variables")

# TODO add ade4

# PCA with pcaMethods::pca
suppressPackageStartupMessages(library("pcaMethods"))
pcaMsvd <- pca(d, method="svd", scale="uv", nPcs=3)
obsMsvd <- fortify(pcaMsvd, type="observations")
varMsvd <- fortify(pcaMsvd, type="variables")


test_that("eigenvalues are equal", {
  
  eigS <- eigenvalues(pcaS)
  eigF <- eigenvalues(pcaF)
  eigM <- eigenvalues(pcaM)

  expect_equivalent(eigS, eigF)
  expect_equivalent(eigS, eigMsvd)
})


test_that("scores are equivalent, proportional to each other", {
  scoresS    <- obsS[,c(".PC1", ".PC2")]
  scoresF    <- obsF[,c(".PC1", ".PC2")]
  scoresMsvd <- obsMsvd[,c(".PC1", ".PC2")]

  # FactoMineR scales scores differently from prcomp, test that the ratio is always the same
  expect_true(length( unique( round( unlist( scoresS / scoresF ), 10) ) ) == 1)
  expect_equivalent(scoresS, scoresMsvd)
})

test_that("loadings are equal or at least proportional to each other", {
  loadingsS    <- varS[,c(".PC1", ".PC2")]
  loadingsF    <- varF[,c(".PC1", ".PC2")]
  loadingsMsvd <- varMsvd[,c(".PC1", ".PC2")]

  expect_equivalent(loadingsS, loadingsF)
  expect_equivalent(loadingsS, loadingsMsvd)
})

test_that("cos2 of observations are equal", {
  expect_equal(obsS$.cos2, obsF$.cos2)
  expect_equal(obsS$.cos2, obsMsvd$.cos2)
})

test_that("cos2 of variables are equal", {
  expect_equal(varS$.cos2, varF$.cos2)
  expect_equal(varS$.cos2, varMsvd$.cos2)
})

test_that("contributions of observations are equal", {
  # FactoMineR scales contributions differently from prcomp, test that the ratio is always the same
  expect_true(length( unique( round( obsS$.contrib / obsF$.contrib, 10) ) ) == 1)
  expect_equal(obsS$.contrib, obsMsvd$.contrib)
})

test_that("contributions of variables are equal or at least proportional to each other", {
  expect_equal(varS$.contrib, varF$.contrib)
  expect_equal(varS$.contrib, varMsvd$.contrib)
})
