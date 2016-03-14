library("testthat")

context("Homogenization of ordination output")

# create some data
set.seed(123)
d <- data.frame(a=runif(5, 1, 10), b=runif(5, 1, 100), c=runif(5, 1, 20))
d <- round(d)

# PCA
pcaS <- stats::prcomp(d, scale=TRUE)
pcaF <- FactoMineR::PCA(d, graph=FALSE, ncp=3)
pcaV <- vegan::rda(d, scale=TRUE)
pcaA <- ade4::dudi.pca(d, scannf=FALSE, nf=3)
pcaM <- pcaMethods::pca(d, method="svd", scale="uv", nPcs=3)

# CA
caF <- FactoMineR::CA(d, ncp=3, graph=FALSE)
caM <- MASS::corresp(d, nf=3)
caC <- ca::ca(d, nd=3)

test_that("eigenvalues are equal", {
  # use stats as reference
  eig_pcaS <- eigenvalues(pcaS)
  expect_equivalent(eig_pcaS, eigenvalues(pcaF))
  expect_equivalent(eig_pcaS, eigenvalues(pcaV))
  expect_equivalent(eig_pcaS, eigenvalues(pcaA))
  expect_equivalent(eig_pcaS, eigenvalues(pcaM))
  
  # use FactoMineR as reference
  eig_caF <- eigenvalues(caF)
  expect_equivalent(eig_caF, eigenvalues(caM))
  expect_equivalent(eig_caF, eigenvalues(caC))
})





# # PCA
# pcaS <- stats::prcomp(d, scale=TRUE)
# pcaF <- FactoMinR::PCA(d, graph=FALSE, ncp=2)
# pcaV <- vegan::rda(d, scale=TRUE)
# pcaA <- ade4::dudi.pca(d, scannf=FALSE, nf=2)
# pcaM <- pcaMethods::pca(d, method="svd", scale="uv", nPcs=2)
#
# # CA
# caF <- FactoMineR::CA(d, ncp=2, graph=FALSE)
# caM <- MASS::corresp(d, nf=3)
# caC <- ca::ca(d, nd=2)
#
#
# test_that("scores are equivalent, proportional to each other", {
#   scoresS    <- obsS[,c(".PC1", ".PC2")]
#   scoresF    <- obsF[,c(".PC1", ".PC2")]
#   scoresMsvd <- obsMsvd[,c(".PC1", ".PC2")]
#
#   # FactoMineR scales scores differently from prcomp, test that the ratio is always the same
#   expect_true(length( unique( round( unlist( scoresS / scoresF ), 10) ) ) == 1)
#   expect_equivalent(scoresS, scoresMsvd)
# })
#
# test_that("loadings are equal or at least proportional to each other", {
#   loadingsS    <- varS[,c(".PC1", ".PC2")]
#   loadingsF    <- varF[,c(".PC1", ".PC2")]
#   loadingsMsvd <- varMsvd[,c(".PC1", ".PC2")]
#
#   expect_equivalent(loadingsS, loadingsF)
#   expect_equivalent(loadingsS, loadingsMsvd)
# })
#
# test_that("cos2 of observations are equal", {
#   expect_equal(obsS$.cos2, obsF$.cos2)
#   expect_equal(obsS$.cos2, obsMsvd$.cos2)
# })
#
# test_that("cos2 of variables are equal", {
#   expect_equal(varS$.cos2, varF$.cos2)
#   expect_equal(varS$.cos2, varMsvd$.cos2)
# })
#
# test_that("contributions of observations are equal", {
#   # FactoMineR scales contributions differently from prcomp, test that the ratio is always the same
#   expect_true(length( unique( round( obsS$.contrib / obsF$.contrib, 10) ) ) == 1)
#   expect_equal(obsS$.contrib, obsMsvd$.contrib)
# })
#
# test_that("contributions of variables are equal or at least proportional to each other", {
#   expect_equal(varS$.contrib, varF$.contrib)
#   expect_equal(varS$.contrib, varMsvd$.contrib)
# })
