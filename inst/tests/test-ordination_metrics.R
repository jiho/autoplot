library("testthat")

context("Homogenization of ordination output")

# create some data
set.seed(123)
n <- 10
d <- data.frame(
  a=runif(n, 1, 2),
  x1=runif(n, 1, 10),
  x2=runif(n, 1, 10),
  x3=runif(n, 1, 20),
  x4=runif(n, 1, 30),
  x5=runif(n, 1, 40),
  x6=runif(n, 1, 50)
)
d <- round(d)
d$a <- letters[d$a]
dnum <- d[-1]

# PCA
# basic and advanced (sup rows, variables, not all components kept, etc)
pcaS <- stats::prcomp(d[,-1], scale=TRUE)
pcaF <- FactoMineR::PCA(d[,-1], graph=FALSE, ncp=6)
pcaFadv<- FactoMineR::PCA(d, graph=FALSE, ncp=2, quali.sup=1, quanti.sup=2, ind.sup=1)
pcaV <- vegan::rda(d[,-1], scale=TRUE)
pcaA <- ade4::dudi.pca(d[,-1], scannf=FALSE, nf=6)
pcaAadv <- ade4::dudi.pca(d[,-1], scannf=FALSE, nf=2)
pcaM <- pcaMethods::pca(d[,-1], method="svd", scale="uv", nPcs=6)
pcaMadv <- pcaMethods::pca(d[,-1], method="svd", scale="uv", nPcs=2)

# CA
# same as above
caF <- FactoMineR::CA(d[,-1], ncp=6, graph=FALSE)
caFadv <- FactoMineR::CA(d[,-1], ncp=2, graph=FALSE, row.sup=1, col.sup=1, quanti.sup=2)
# caFadv <- FactoMineR::CA(d, ncp=2, graph=FALSE, row.sup=1, quali.sup=1, col.sup=2, quanti.sup=3)
# NB: bug in FactoMineR
caM <- MASS::corresp(d[,-1], nf=6)
caMadv <- MASS::corresp(d[,-1], nf=2)
caC <- ca::ca(d[,-1], nd=6)
caCadv <- ca::ca(d[,-1], nd=2, suprow=2, supcol=1)

test_that("eigenvalues are equal among packages", {
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

test_that("number of active data rows is correctly computed", {
  n <- nrow(d)
  
  expect_equal(n, nr(pcaS))
  expect_equal(n, nr(pcaF))
  expect_equal(n, nr(pcaV))
  expect_equal(n, nr(pcaA))
  expect_equal(n, nr(pcaM))

  expect_equal(n, nr(caF))
  expect_equal(n, nr(caM))
  expect_equal(n, nr(caC))

  # advanced cases
  expect_equal(n-1, nr(pcaFadv))  # one supplementary observation
  expect_equal(n, nr(pcaAadv))
  expect_equal(n, nr(pcaMadv))

  expect_equal(n-1, nr(caFadv))   # idem
  expect_equal(n, nr(caMadv))
  expect_equal(n-1, nr(caCadv))   # idem
})

test_that("number of active columns is correctly computed", {
  n <- ncol(d)
  
  # basic cases
  expect_equal(n-1, nc(pcaS))    # 1 non-numeric column
  expect_equal(n-1, nc(pcaF))    # idem
  expect_equal(n-1, nc(pcaV))    # idem
  expect_equal(n-1, nc(pcaA))    # idem
  expect_equal(n-1, nc(pcaM))    # idem

  expect_equal(n-2, nc(caF))     # 1 non-numeric column, 1-less dim for CA
  expect_equal(n-2, nc(caM))     # idem
  expect_equal(n-2, nc(caC))     # idem

  # advanced cases
  expect_equal(n-2, nc(pcaFadv)) # 1 quali.sup, 1 quanti.sup
  expect_equal(n-1, nc(pcaAadv)) # 1 non-numeric column
  expect_equal(n-1, nc(pcaMadv)) # idem

  expect_equal(n-4, nc(caFadv))  # 1 non-numeric, 1 col.sup, 1 quanti.sup, 1-less dim for CA
  expect_equal(n-2, nc(caMadv))  # 1 non-numeric, 1-less dim for CA
  expect_equal(n-3, nc(caCadv))  # 1 non-numeric, 1 supcol, 1-less dim for CA
})

test_that("number of dimensions kept is correctly computed", {
  n <- ncol(d)
  
  # basic cases
  expect_equal(n-1, npc(pcaS))   # 1 non-numeric column
  expect_equal(n-1, npc(pcaF))   # idem
  expect_equal(n-1, npc(pcaV))   # idem
  expect_equal(n-1, npc(pcaA))   # idem
  expect_equal(n-1, npc(pcaM))   # idem
                                 
  expect_equal(n-2, npc(caF))    # 1 non-numeric column, 1-meaningless dim for CA
  expect_equal(n-2, npc(caM))    # idem (NB: we explicitly exclude it for MASS)
  expect_equal(n-2, npc(caC))    # idem

  # advanced cases
  expect_equal(2, npc(pcaFadv))  # only 2 kept
  expect_equal(2, npc(pcaAadv))  # idem
  expect_equal(2, npc(pcaMadv))  # idem
                                 
  expect_equal(2, npc(caFadv))   # idem
  expect_equal(2, npc(caMadv))   # idem
  expect_equal(2, npc(caCadv))   # idem
})

test_that("scale 0 scores are equal among packages", {
  abs_scores <- function(x, ...) {
    sco <- scores(x, scaling=0, ...)
    abs(sco[sapply(sco, is.numeric)])
  }
  
  # use stats as reference
  scores_pcaS <- abs_scores(pcaS, which="row")
  expect_equivalent(scores_pcaS, abs_scores(pcaF, which="row"))
  expect_equivalent(scores_pcaS, abs_scores(pcaV, which="row"))
  expect_equivalent(scores_pcaS, abs_scores(pcaA, which="row"))
  expect_equivalent(scores_pcaS, abs_scores(pcaM, which="row"))
  scores_pcaS <- abs_scores(pcaS, which="col")
  expect_equivalent(scores_pcaS, abs_scores(pcaF, which="col"))
  expect_equivalent(scores_pcaS, abs_scores(pcaV, which="col"))
  expect_equivalent(scores_pcaS, abs_scores(pcaA, which="col"))
  expect_equivalent(scores_pcaS, abs_scores(pcaM, which="col"))
  
  # use FactoMineR as reference
  scores_caF <- abs_scores(caF, which="row")
  expect_equivalent(scores_caF, abs_scores(caM, which="row"))
  expect_equivalent(scores_caF, abs_scores(caC, which="row"))
  scores_caF <- abs_scores(caF, which="col")
  expect_equivalent(scores_caF, abs_scores(caM, which="col"))
  expect_equivalent(scores_caF, abs_scores(caC, which="col"))
})

test_that("scaling matches FactoMineR", {
  rows <- data.frame(abs(pcaF$ind$coord))
  cols <- data.frame(abs(pcaF$var$coord))
  
  abs_scores <- function(x, ...) {
    sco <- scores(x, ...)
    abs(sco[sapply(sco, is.numeric)])
  }
  
  for (pca_obj in list(pcaS, pcaF, pcaV, pcaA, pcaM)) {
    expect_equivalent(rows, abs_scores(pca_obj, which="rows", scaling=1)[,1:6])
    expect_equivalent(cols, abs_scores(pca_obj, which="columns", scaling=2)[,1:6])    
  }
})


test_that("cos2 matches FactoMineR", {
  rows <- pcaF$ind$cos2
  for (pca_obj in list(pcaS, pcaF, pcaV, pcaA, pcaM)) {
    for (i in 1:6) {
      expect_equivalent(rows[,i], augment(pca_obj, which="row", dim=i)$.cos2)
    }
  }
  cols <- pcaF$var$cos2
  for (pca_obj in list(pcaS, pcaF, pcaV, pcaA, pcaM)) {
    for (i in 1:6) {
      expect_equivalent(cols[,i], augment(pca_obj, which="col", dim=i)$.cos2)
    }
  }
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
