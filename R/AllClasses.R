## seed set class
#' \code{seedset-class}
#'
#' Object representing the seed sets of a given metabolic network
#'
#' \section{Slots}{
  #' \describe{
  #'    \item{seeds}{a character list represents seeds of a given metabolic network which
  #'                is composed of  the KEGG compound index.}
  #'    \item{confidenceScore}{a numeric list represents confidence sore of the seed}
  #'  }
#' @name seedset-class
#' @rdname seedset-class
#' @exportClass seedset
#'

### Validity check
.check.seedset<- function(object)
{
  if(!is(object, "seedset")) stop("object has to be of class \"seedset\" ")

  errors <- character()

  if(!is.list(object@seeds))
    errors <- c(errors, "seeds must be a list")
  if(!is.list(object@confidenceScore))
    errors <- c(errors, "confidenceScore must be a list")
  if (length(object@seeds) != length((object@confidenceScore)))
    errors <- c(errors, "length of seeds and confidenceScore should be equal ")
  if(length(errors) == 0){
    TRUE
  }else{
    errors
  }
}

## set class
setClass("seedset", slot=list(seeds="list",confidenceScore="list"),
         prototype = list(seeds=NULL, confidenceScore=NULL),
         validity = .check.seedset, sealed=TRUE)
