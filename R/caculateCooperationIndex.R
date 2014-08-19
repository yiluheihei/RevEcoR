##==============================================================================
##==============================================================================
#' Caculating the metabolic complementarity index
#' 
#' Caculating the metabolic complementarity index of g1 in the presence of g2
#' 
#'@param g1, igraph object, a metabolic network
#'
#'@param g2, igraph object, a metabolic network, the complementary network of g1
#'
#'@param threshold, the cutoff of confidence score to be serve as a seed set, 
#'default is 0
#'
#'@return a number, range from 0 to 1

complementarityIndex <- function(g1,g2, threshold=0){
  if (!is.igraph(g1) || !is.igraph(g2))
    stop("Both g1 and g2 must be igraph object")
  seed.set1 <- getSeedSets(g1,threshold)@seeds
  seed.set2 <- getSeedSets(g2, threshold)@seeds
  nonseed2 <- setdiff(V(g2)$name, unlist(seed.set2))
  complement.cpd <- setdiff(unlist(seed.set1),unlist(seed.set2)) %>%
    intersect(., nonseed2)
  if(length(complement.cpd)){
    norm.seed <- lapply(complement.cpd,function(x)lapply(seed.set1,
      function(y)match(x,y,nomatch=0))) %>%
      lapply(.,function(x)which(x>0)) %>%
      unique %>%
      length
    complementary.index <- norm.seed/length(seed.set1)
  }else{
    complementary.index <- 0
  }
    return(complementary.index)
}

##==============================================================================
##==============================================================================

#' Caculating the metabolic competition index
#'  
#' Caculating the metabolic competition index of g1 in the presence of g2
#' 
#'@param g1, igraph object, a metabolic network
#'
#'@param g2, igraph object, a metabolic network, the complementary network of g1
#'
#'@param threshold, the cutoff of confidence score to be serve as a seed set, 
#'default is 0.2


competitionIndex <- function(g1,g2, threshold=0){
  if (!is.igraph(g1) || !is.igraph(g2))
    stop("Both g1 and g2 must be igraph object")
  seed.set1 <- getSeedSets(g1,threshold)@seeds
  seed.set2 <- getSeedSets(g2, threshold)@seeds
  intersect.seed <- intersect(unlist(seed.set1),unlist(seed.set2))
  if(length(intersect.seed)){
    norm.seed <- lapply(intersect.seed,function(x)lapply(seed.set1,
      function(y)match(x,y,nomatch=0))) %>%
      lapply(.,function(x)which(x>0)) %>%
      unique %>%
      length
    competition.index <- norm.seed/length(seed.set1)
  }else{
    competition.index <- 0
  }
  return(competition.index)
}

##==============================================================================
##==============================================================================
BSIscore <- function(g1, g2, threshold=0){
  if (!is.igraph(g1) || !is.igraph(g2))
    stop("Both g1 and g2 must be igraph object")
  seed.set1 <- getSeedSets(g1,threshold)@seeds
  seed.set2 <- getSeedSets(g2, threshold)@seeds
  intersect.seed <- intersect(unlist(seed.set1), V(g2)$name)
  if(length(intersect.seed)){
    norm.seed <- lapply(intersect.seed,function(x)lapply(seed.set1,
      function(y)match(x,y,nomatch=0))) %>%
      lapply(.,function(x)which(x>0)) %>%
      unique %>%
      length
    BSI.score <- norm.seed/length(seed.set1)
  }else{
    BSI.score <- 0
  }
  return(BSI.score) 
}

##==============================================================================
##==============================================================================
#'Caculating the metabolic competition and complementarity index
#'
#'Caculating the metabolic competition complementarity index among all metabolic 
#'networks
#'
#'@param g, igraph that represents a metabolic network, see \code{\link{reconstructGsMN}}
#'
#'@param ..., a list of metabolic networks or a network append to g
#'
#'@param threshold threshold, the cutoff of confidence score to be serve as a 
#'seed set, default is 0.2
#'
#'@export
#'
#'@details Metabolic competition index is defined as the fraction of compounds 
#'in a species seed set of metabolic network that are alse included in its 
#'partner; However, metabolic complementarity index is the fraction of 
#'compounds in one species seed set of metabolic network appearing in the 
#'metabolic network but not in the seed set of its partner; The biosynthetic 
#'support score represents the extent to which the metabolic requirements of a 
#'potential parasitic organism can be supported by the biosynthetic capacity of
#'a potential host. It is measured by calculating the fraction of the source 
#'components of a, in which at least one of the compounds can be found in the 
#'network of b. However, seed compounds are associated with a confidence score 
#'(1/size of SCC), so this fraction is calculated as a mormalized weighted sum. 
#'
#'The ith row and jth col elements of the returnd matrix represents the 
#'metabolic competition index or complementarity index of the ith network on the
#'jth metabolic network.
#'
#'@return a cooperation index matrix whose nrow and ncol is equal to the number 
#'of species to be compared, for more see details.

caculateCooperationIndex <- function(g, ...,threshold=0){
  g <- list.append(g, ...)
  if (length(g)<2)
    stop("At least two species to compare")
  competition.index <- matrix(0,length(g),length(g))
  complementarity.index <- matrix(0,length(g),length(g))
  bsi.score <- matrix(0,length(g),length(g))
  index <- permutations(length(g),2)
  competition.v<- apply(index,1,function(x)g[x]) %>%
    sapply(., function(x)competitionIndex(x[[1]],x[[2]],threshold))
  for(i in 1:nrow(index))
    competition.index[index[i,1],index[i,2]] = competition.v[i]
  complementarity.v<- apply(index,1,function(x)g[x]) %>%
    sapply(., function(x)complementarityIndex(x[[1]],x[[2]],threshold))  
  for(i in 1:nrow(index))
    complementarity.index[index[i,1],index[i,2]] = complementarity.v[i]
  bsi.v<- apply(index,1,function(x)g[x]) %>%
    sapply(., function(x)BSIscore(x[[1]],x[[2]],threshold))  
  for(i in 1:nrow(index))
  bsi.score[index[i,1],index[i,2]] = bsi.v[i]
  #row.names(interaction.index) <- paste0(rep("g",length(g)),seq(length(g)))
  #colnames(interaction.index) <- row.names(interaction.index)
  return(list(competition.index = competition.index, 
    complementarity.index = complementarity.index, bsi.score = bsi.score))
}
##==============================================================================
##==============================================================================