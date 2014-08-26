#'@title Reconstuction of the specific-organism genome-scale metabolic network
#'
#'@description Reconstruction of genome-scale metabolic network (GsMN) whose
#'nodes represents compounds and whose edges represents reactions.
#'
#'@param metabolic.data, df, more details see function \code{getOrgMetabolicData}
#'and \code{details}
#'
#'@details The input of this function can be of two forms. If organims is 
#'collected in KEGG database, it can be obtained with \code{getOrgMetabolicData}. 
#'If not, \code{metabolic.data} could be a  table format file which contains the
#'KEGG Orthology annotated information on this organism, e.g. we can download 
#'this KO annotation profile in the \url{https://img.jgi.doe.gov} website for 
#'species detected in a human microbime which not contained in KEGG organism 
#'database. Several functions, such as \code{link{read.table}} and 
#'\code{\link{read.delim}} could help us to read KO annotation profile.
#'
#'@return igraph object
#'
#'@export
#'
#'@seealso \code{\link{getOrgMetabolicData}}
#'
#'@examples
#'## not run (organism in KEGG)
#'## metabolic.data <- getOrgMetabolicData("buc", refresh = TRUE)
#'## g <- reconstructGsMN(metabolic.data)
#'
#'## species detected in a human microbiome
#' annodir <- system.file("extdata","koanno.tab",package = "ReverseEcologyR")
#' metabolic.data <- read.delim2(file=annodir,stringsAsFactors=FALSE)
#' g2 <- reconstructGsMN(metabolic.data)

reconstructGsMN <- function(metabolic.data){
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
  if (match(".attrs.name",names(metabolic.data),nomatch=0)){
    message("metabolic.data is the reaction information from KEGG organism pathway map...")
    metabolites <- metabolic.data[,c(2,3)]
  }else{
    message("metabolic.data is the KEGG Orthology annotation profile of the species...")
    if (FALSE){
      if (!file.exists(file.path(Sys.getenv("home"),".mmnet","RefDbcache.rda"))){
        message("Update the KEGG pathway data and save local...")
        updateKEGGPathway()
      }
      load(file.path(Sys.getenv("home"),".mmnet","RefDbcache.rda"))
    }
    data(RefDbcache,package="mmnet")
    ko <- intersect(RefDbcache$ko, metabolic.data[[1]])
    index <- match(ko,RefDbcache$ko)
    metabolites <- RefDbcache[c("substrate","product")] %>%
      lapply(.,"[",index)
    metabolites <- data.frame(substrate=do.call(cbind,metabolites[1]),
      product=do.call(cbind,metabolites[2]),row.names=NULL,stringsAsFactors=FALSE)
  }
  node <- unique(unlist(metabolites))
  net.matrix <- matrix(0,length(node),length(node))
  ##test = lapply(node,function(x)lapply(metabolites[,1],intersect,x))
  row.index <- lapply(metabolites[,1],match,node)
  col.index <- lapply(metabolites[,2],match,node)
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




