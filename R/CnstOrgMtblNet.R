#'@title reconstuction of the specific-organism metabolic network
#'
#'@param pathway.info, df, more details see function "GetOrgPathwayInfo.R"
#'
#'@return igraph object
#'
#'@export

CnstOrgMtblNet <- function(pathway.info){
  metabolites <- pathway.info[,c(2,3)]
  node <- unique(unlist(metabolites))
  net.matrix <- matrix(0,length(node),length(node))
  ##test = lapply(node,function(x)lapply(metabolites[,1],intersect,x))
  row.index <- lapply(pathway.info[,2],match,node)
  col.index <- lapply(pathway.info[,3],match,node)
  row.col <- mapply(expand.grid,row.index,col.index)
  row.index <- unlist(row.col[1,])
  col.index <- unlist(row.col[2,])
  #mapply(function(x,y)net.matrix[x,y]=1,row.index,col.index)
  for (i in 1:length(row.index))
    net.matrix[row.index[i],col.index[i]] <- 1
  diag(net.matrix) <- 0
  g <- graph.adjacency(net.matrix)
  V(g)$name <- node
  ## omit the glycans
  node2 <- str_count(node, "^gl") %>%
    is_greater_than(0) %>%
    extract(node,.)
  ## drop the small dis-connetect components
  g <- delete.vertices(g, node2)
  #g.community <- clusters(g)
  #community.size <- sizes(g.community) 
  #small.size <- which(community.size <= 10)
  #small.vertex <- match(membership(g.community), small.size)
  #small.vertex <- V(g)$name[!is.na(small.vertex)]  
  #g <- delete.vertices(g, small.vertex)
}
