# Get number of potential objects = *active* data rows
# generic
nr <- function(x) { UseMethod("nr") }
# define methods
nr.prcomp <- function(x) { nrow(x$x) }
nr.PCA    <- function(x) { nrow(x$ind$coord) }
nr.rda    <- function(x) { nrow(x$CA$u) }
nr.pca    <- function(x) { nrow(x$li) }
nr.pcaRes <- function(x) { x@nObs }

nr.CA <- function(x) { nrow(x$row$coord) }
nr.correspondence <- function(x) { nrow(x$rscore) }
nr.ca <- function(x) { sum(!is.na(x$rowinertia)) }


# Get number of potential dimensions = *active* data columns
# generic
nc <- function(x) { UseMethod("nc") }
# define methods
nc.prcomp <- function(x) { length(x$sdev) }
nc.PCA    <- function(x) { nrow(x$eig) }
nc.rda    <- function(x) { length(x$CA$eig) }
nc.pca    <- function(x) { x$rank }
nc.pcaRes <- function(x) { x@nVar }

# NB for CA, the number of potential dimensions is ncol - 1
nc.CA <- function(x) { nrow(x$eig) }
nc.correspondence <- function(x) { ncol(x$Freq) - 1 }
nc.ca <- function(x) { sum(!is.na(x$colinertia)) - 1 }


# Get number of kept principal components
# generic
npc <- function(x) { UseMethod("npc") }
# define methods
npc.prcomp <- nc.prcomp
npc.PCA    <- function(x) { x$call$ncp }
npc.rda    <- nc.rda
npc.pca    <- function(x) { x$nf }
npc.pcaRes <- function(x) { x@nPcs }

npc.CA <- function(x) { x$call$ncp }
npc.correspondence <- function(x) { length(x$cor) }
npc.ca <- function(x) { length(x$sv) }
