get.targets <- function(x, Names, pval=0.001){
	bound.genes <- Names[(x < pval) & !is.na(x)]
	return(as.character(bound.genes))
} 

get.all.targets <- function(tfs, Names, pval=0.001){
	bound.targets <- list()
	tf.names <- names(tfs)
	
	for(name in tf.names){
		targets <- get.targets(tfs[[name]], Names, pval)
		bound.targets[[name]] <- targets
	}
	return(bound.targets)
}

targets.in.common <- function(g, name1, name2){
	targets1 <- neighbors(g, name1)
	targets2 <- neighbors(g, name2)
	common.targets <- intersect(targets1, targets2)
	return(common.targets)
}

number.in.common <- function(g, tfs){
  # not very efficient, but easy to understand
  n <- length(tfs)
  m <- mat.or.vec(n, n)
  rownames(m) <- tfs
  colnames(m) <- tfs
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      m[i,j] <- length(targets.in.common(g, tfs[i], tfs[j]))
      m[j,i] <- m[i,j]
    }
  }
  return(m)
}

giant.component <- function(g, ...) {
  cl <- clusters(g, ...)
  induced_subgraph(g,  which(cl$membership == which.max(cl$csize)))
}
