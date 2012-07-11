context("Get data from a PCA object")

library("FactoMineR")
library("pcaMethods")

# create some data
set.seed(123)
d <- data.frame(a=runif(5), b=runif(5), c=runif(5))


test_that("fortify does concatenate data", {
  pca <- prcomp(d)
  expect_equal(fortify(pca, data=d)[1:ncol(d)], d)
})

test_that("fortify concatenates all data, including columns not used in the PCA", {
  pca <- prcomp(d[1:2])
  expect_equal(fortify(pca, data=d)[1:ncol(d)], d)
})

test_that("fortify concatenates all data, including rows not used in the PCA", {
  pca <- prcomp(d[1:4,])
  expect_equal(fortify(pca, data=d)[1:ncol(d)], d)
  expect_true(is.na(fortify(pca, data=d)$.PC1[5]))

  pca <- PCA(d[1:4,], graph=F)
  expect_equal(fortify(pca, data=d)[1:ncol(d)], d)
  expect_true(is.na(fortify(pca, data=d)$.PC1[5]))

  pca <- pca(d[1:4,], graph=F)
  expect_equal(fortify(pca, data=d)[1:ncol(d)], d)
  expect_true(is.na(fortify(pca, data=d)$.PC1[5]))
})

test_that("row order does not affect fortify", {
  rd <- d[sample.int(nrow(d)),]

  pca <- prcomp(d)
  rpca <- prcomp(rd)
  expect_equal(abs(fortify(pca, data=d)$.PC1), abs(fortify(rpca, data=d)$.PC1))

  pca <- PCA(d, graph=FALSE)
  rpca <- PCA(rd, graph=FALSE)
  expect_equal(abs(fortify(pca, data=d)$.PC1), abs(fortify(rpca, data=d)$.PC1))

  pca <- pca(d)
  rpca <- pca(rd)
  expect_equal(abs(fortify(pca, data=d)$.PC1), abs(fortify(rpca, data=d)$.PC1))
})

# TODO test using imputed columns rather than data columns
