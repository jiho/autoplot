# Extract eigenvalue diagnostics from an ordination analysis
tidy_ordination <- function(x, ...) {
  # extract eigenvalues
  eig <- eigenvalues(x)
  
  # test that all potential PCs are there to compute proportions
  npc_max <- nc(x)
  if (length(eig) == npc_max) {
    prop.variance <- eig/sum(eig)
    cum.prop.variance <- cumsum(prop.variance)
  } else {
    prop.variance <- NA
    cum.prop.variance <- NA
  }
  
  # format result as a data.frame
  data.frame(
    term=paste0("PC", 1:length(eig)),
    eigenvalues=eig,
    prop.variance,
    cum.prop.variance
  )
}
