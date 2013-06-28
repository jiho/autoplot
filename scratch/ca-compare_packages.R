source("fortify-fact.R")
source("autoplot-fact.R")

# FactoMineR is the reference since it seems to be the most capable and best documented package of the lot

## CA                                                       {

# with MASS
caM = corresp(d, nf=2)

# with FactoMineR
caF = CA(d)

# with ade4

# Checks

# eigenvalues
(eig = caM$cor^2)
caF$eig$eigenvalue

# column scores
(cscores = t(t(caM$cscore) * caM$cor))
caF$col$coord

# row scores
(rscores = t(t(caM$rscore) * caM$cor))
caF$row$coord

# cos2
cscores^2 / rowSums(cscores^2)
caF$col$cos2

rscores^2 / rowSums(rscores^2)
caF$row$cos2

# contribution
t(t(cscores^2*(colSums(d)/sum(d)))/eig) * 100
caF$col$contrib

t(t(rscores^2*(rowSums(d)/sum(d)))/eig) * 100
caF$row$contrib

fortify(caM, data=d)
fortify(caF)

autoplot(caM, data=d)
autoplot(caF)
autoplot(caM, data=d, mapping=aes(size=.contrib, alpha=.cos2))
autoplot(caF, mapping=aes(size=.contrib, alpha=.cos2))

# }


## MCA                                                      {

# }


