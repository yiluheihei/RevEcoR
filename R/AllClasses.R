.check.seedset<- function(object)
{
  if(!is(object, "seedset")) stop("object has to be of class \"seedset\" ")
  
  errors <- character()
  if (!is.igraph(object@GsMN))
    errors <- c(errors, "GsMN must be a igrpah object")  
  if(!is.list(object@seeds))
    errors <- c(errors, "seeds must be a list")
  if(length(errors) == 0){
    TRUE
  }else{
    errors
  }
}
#' \code{seedset-class}
#'
#' Object representing the seed sets of a given metabolic network
#'
#'@slot GsMN, a igraph network
#'@slot seeds, a character list represents seeds of a given metabolic 
#'      network which is composed of  the KEGG compound index.
#'@slot confidenceScore, a numeric list represents confidence sore of the 
#'      seed
#'@section method:
#'    \itemize{
#'      \item{getGsMN, \code{signature(object = "seedset")}:
#'        get the genome scale metabolic network whose seed set is caculated}
#'      \item{length, \code{signature(object = "seedset")}:
#'        return the number of source SCC}
#'      \item{seedSize, \code{signature(object = "seedset")}:
#'        returns the sizes of each source SCCs}
#'      \item{nonseed, \code{signature(object = "seedset")}:
#'        the non seeds of the GsMN}
#'      \item{show, \code{signature(object = "seedset")}:
#'        show the short summary of a seedset class}
#'      \item{confidencescore, \code{signature(object = "seedset")}:s
#'        confidence score of the seed set}
#'  }
#' @name seedset-class
#' @rdname seedset-class
#' @exportClass seedset
#' @seealso \code{\link{getSeedSets}}
#' 
setOldClass("igraph")
setClass("seedset", slot=list(GsMN="igraph",seeds="list"),
         prototype = list(GsMN=NULL,seeds=NULL),
         validity = .check.seedset)


#' Size of the each seed source component
#'
#' Caculate the size of each seed source component.
#'
#'@export 
#'@rdname seedSize-methods
#'@param obj, \code{seedset} class
#'@name seedSize-methods
#'@aliases seedSize seedSize-methods 
#'@docType methods
#'@seealso \code{\link{seedset-class}}
#'@return a vector represents size of each source seed componet of network


setGeneric("seedSize",
  function(object){standardGeneric("seedSize")})

setMethod("seedSize",signature="seedset",
  function(object){
    listLen(object@seeds)
  }
)

#' The genome scale metabolic network 
#'
#' T he genome scale metabolic network (GsMN) whose seed set is caculated.
#'
#'@export 
#'@rdname getGsMN-methods
#'@param obj, \code{seedset} class
#'@name getGsMN-methods
#'@aliases getGsMN getGsMN-methods 
#'@docType methods
#'@seealso \code{\link{seedset-class}}
#'@return a igraph


setGeneric("getGsMN",
  function(object, field){standardGeneric("getGsMN")})

setMethod("getGsMN",signature="seedset",
  function(object){
    object@GsMN
  }
)

#' Non seed of the network 
#'
#' Non seed of the network.
#'
#'@export 
#'@rdname nonseed-methods
#'@param obj, \code{seedset} class
#'@name nonseed-methods
#'@aliases nonseed nonseed-methods 
#'@docType methods
#'@seealso \code{\link{seedset-class}}
#'@return a vector


setGeneric("nonseed",
  function(object, field){standardGeneric("nonseed")})

setMethod("nonseed",signature=c(obj="seedset"),
  function(object){
    non.seed  <- V(object@GsMN)$name
    non.seed  <- setdiff(non.seed,unlist(object@seeds))
  }
)


#' Conficence score 
#'
#' Caculate confidence score of seed set 
#'
#'@export 
#'@rdname confidencescore-methods
#'@param obj, \code{seedset} class
#'@name confidencescore-methods
#'@aliases confidencescore confidencescore-methods 
#'@docType methods
#'@seealso \code{\link{seedset-class}}
#'@return a vector


setGeneric("confidencescore",
  function(object, field){standardGeneric("confidencescore")})

setMethod("confidencescore",signature=c(obj="seedset"),
  function(object){
    confidence.score  <- 1/seedSize(object) 
      mapply(function(x,y)rep(x,y),confidence.score,seedSize(object))
  }
)
#'  the length generic function
#'
#' Caculate the number of the seed source components.
#'
#'@export 
#'@rdname length-methods
#'@param obj, \code{length} class
#'@name length-methods
#'@aliases length length-methods 
#'@docType methods
#'@seealso \code{\link{seedset-class}}
#'@return a interger

setMethod("length",signature="seedset",
  function(x){
    length(x@seeds)
  }
)

#' The show generic function
#'
#' Show a short summary of seedset object
#'
#'@export 
#'@rdname showmethods
#'@param obj, \code{seedset} class
#'@name show-methods
#'@aliases show show-methods 
#'@docType methods
#'@seealso \code{\link{seedset-class}}


setMethod("show",signature=c(obj="seedset"),
  function(object){
    cat("Object of class ", class(object), "\n", sep = "")
    igraph::print.igraph(object@GsMN)
    cat("seedset length ", length(object), "\n")
  }
)

