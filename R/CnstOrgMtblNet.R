#'@title reconstuction of the specific-organism metabolic network
#'
#'@param pathway.info, df, more details see function "GetOrgPathwayInfo.R"
#'
#'@return igraph object
#'
#'@export

CnstOrgMtblNet <- function(pathway.info){
  # delete the small group nodes
  deleteSG <- function(g,threshold = 10){ 
  if (!is.igraph(g))
      stop("Not a igraph object")
    delet.node <- NULL
    Step.count <- 6 
    while (all((neighborhood.size(g,Step.count)-neighborhood.size(g,Step.count-1))!=0)){
      Step.count <- Step.count + 1
    }
    delet.node <- which(neighborhood.size(g,Step.count) <= threshold)
    gf <- delete.vertices(g,delet.node)
    return(gf)
  }
  
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
  ## omit the glycans and drugs
  node2 <- str_count(node, "^gl|^dr") %>%
    is_greater_than(0) %>%
    extract(node,.)
  ## drop the small dis-connetect components
  g <- delete.vertices(g, node2) %>%
    deleteSG(threshold = 10)
  g
}
