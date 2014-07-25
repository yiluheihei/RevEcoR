#' Identify seed compounds of each organism
#' 
#' Detect a given metabolic network and idendity the seed compounds of each organism
#' 
#' @param g,  an igraph object which represents a given organism-specific metaboliic network
#' @param threshold, numeric constant ranges from 0 to 1, default is 0.
#' @details All the compound in the same source SCC all equally to be included in the seed set,
#' each of these compounds was assigned a confidence level, C=1/(size of souce SCC), denoting
#' the compounds probability of being a seed. This threshold was used to determin whether a 
#' compound should be a seed.
#' @return seed set compounds of the given organism-specific metabolic network,
#' genome scale metabolic network is also returned
#' @export
#' @seealso \code{\link{KosarajuSCC}}
getSeedSets <- function(g, threshold = 0){
  if (!is.igraph(g))
    stop("Not a igraph object")
  ## check scc is a source seed sets or not 
  checkSCC <- function(g, x){
    edge.in <- lapply(x, subcomponent, g=g, mode="in")
    edge.out <- lapply(x,subcomponent,g=g, mode="out")
    diff.in <- setdiff(unlist(edge.in), unlist(x))
    diff.out <- setdiff(unlist(edge.out), unlist(x))
    if (length(diff.in)==0 && length(diff.out)>0){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  ## get the seed sets of a given network  
  g.scc <- KosarajuSCC(g)
  minsize.scc <- floor(1/threshold)
  index <- which(listLen(g.scc)<=minsize.scc) 
  min.scc <- g.scc[index]
  seeds <- sapply(min.scc,checkSCC,g=g) %>%
    extract(min.scc,.) %>%
    lapply(.,function(x)extract(V(g)$name,x))
  return(list(GsMN=g,seedset=seeds))
}