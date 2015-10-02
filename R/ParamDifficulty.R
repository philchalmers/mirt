# Extract parameter d for data
#   considering collapsing

# d <- matrix[i.item, i.cat]
# R <- matrix[i.subj, i.item]

#d.new <- matrix(Inf, n.item, ncat-1)
#d.new <- cbind(-Inf, d.new)

paramDiff=function(d, R) {

  d.new <- list()

  stopifnot(is.matrix(d) | is.list(d))

  if (is.matrix(d)) {
    d.aug <- cbind(-Inf, d)

    for (i.item in 1:n.item) {
        x <- R[, i.item]
        s <- sort(unique(x))
        cat.exist <- 1:ncat %in% s
        d.new[[i.item]] <- d.aug[i.item, cat.exist][-1]
        names(d.new[[i.item]])=paste("d",1:length(d.new[[i.item]]), sep="")
    }}
  if (is.list(d)) {# Yet to be checked!
      for (i.item in 1:n.item) {
          x <- R[, i.item]
          s <- sort(unique(x))
          cat.exist <- 1:ncat %in% s
          d.aug <- c(-Inf, d[[i.item]])
          d.new[[i.item]] <- d.aug[cat.exist][-1]
          names(d.new[[i.item]])=paste("d",1:length(d.new[[i.item]]), sep="")
      }}
  return(d.new)
}









