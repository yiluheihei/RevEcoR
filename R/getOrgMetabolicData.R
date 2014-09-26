#' Get organism metabolic data from KEGG database
#' 
#' This function helps us  to obtain the specific-organism pathway map, prasing 
#' this maps to get metabolic data contains reaction, substrate and product.
#' 
#' @param org, three characters, the KEGG organism code, e.g. "buc".
#' 
#' @param path, character which specify the dir for saving the local metabolic 
#' data, default is your home dir of R. More information see \code{details}.
#' 
#' @param refresh, logical, whether refresh the exsits metabolic data from KEGG,
#' detault is FALSE.
#' 
#' @details \code{getOrgMetabolicData} helps us to obtain metabolic data of a 
#' given organism from KEGG database. However, keeping a downloaded metabolic 
#' data in local is necessary. It not only can help keeping the data accurate 
#' and speed up the process of network reconstruction, but also allows you to 
#' reuse the data for other purposes.
#' 
#' Futher more, the reference data such as the KO reference hierarchy does not 
#' change frequently.Neither does themetabolic data of organisms whose genome is
#' sequenced very early. However, to those whose genome is undersequenced or just
#' got completely annotated completely not very long ago, their metabolic data 
#' is updated frequently. In this case, assigning different TTL (Time-to-live) to
#' different types of data will reduce network overhead and data retrieval
#' time caused by remotely visiting KEGG frequently. For example, it would be 
#' helpful to set the renewal period of the KO reference hierarchy data as 30 
#' days while that of the organism-specific.
#' 
#' As \code{getOrgMetabolicData}, check for the existence of queried orgaism 
#' metabolic data in the local database. If the data is not there, retrieve it 
#' from KEGG and write into the local database. If not, return the data 
#' immediately from the local database.
#' 
#' @export
#' 
#' @return a three length df, consists of  reaction name, substrates and products

getOrgMetabolicData <- function(org, path = Sys.getenv("home"), refresh = FALSE){
  ## list of orgnism-specific kegg pathway map 
  org.pathway <- keggList("pathway",org) %>%
    names(.) %>%
    sapply(.,str_replace_all,"path:","")

  ##-----------------------------------------------------------##

  ## information of reaction,direction,sub-pro for pathway map ##
  PathwayInfo <- function(pathway){
    kgml <- keggGet(pathway,"kgml")
    reaction <- lapply(getNodeSet(xmlParse(kgml), "//reaction"),xmlToList)
    if(length(reaction)){
      reaction <- lapply(reaction,unlist)
      attr <- c(".attrs.name","substrate.name","product.name",".attrs.type")
      reaction.info <- lapply(attr,function(x)lapply(reaction,
        function(y)y[names(y)==x]))
      return(reaction.info)
    }else{
      return(list())
    }
  }
  
  ##download the metabolic data from KEGG 
  refreshData <- function(org.pathway){
    ##------integrate all pathway information for each specie------##
    pathway.info <- lapply(org.pathway,PathwayInfo)
    pathway.info <- pathway.info[listLen(pathway.info)>0]
    ##-------------------multiple rn-------------------##
    ##rn.name =sapply(rn.name,str_replace_all,"rn:","")
    rn.info <- lapply(c(1,2,3,4),function(x)Reduce("c",
      lapply(pathway.info,"[[",x))) %>%
      do.call(rbind, .) %>%
      t
    rn.name <- rn.info[,1]
    multi.index <- which(sapply(rn.name,nchar)>9)
    multi.rn <- sapply(rn.name[multi.index],strsplit,"\\s",perl=TRUE)
    multi.len <- listLen(multi.rn)
    multi.metabolites <- rn.info[multi.index,-1]
    multi.metabolites <- multi.metabolites[rep(1:nrow(multi.metabolites),
      times=multi.len),]
    multi.info <- cbind(unlist(multi.rn),multi.metabolites)
    row.names(multi.info) <- unlist(multi.rn)
    rn.info2 <- rbind(rn.info[-multi.index,],multi.info)
    colnames(rn.info2) <- c(".attrs.name","substrate.name","product.name",".attrs.type")
    ##-------------------------------------------------##
    ## nodes represents multiple cpds
    rn.info2[,2] <- rn.info2[,2] %>%
      sapply(strsplit,"\\s") %>%
      lapply(unlist)
    rn.info2[,3] <- rn.info2[,3] %>%
      sapply(strsplit,"\\s") %>%
      lapply(unlist) 
    
    ##----------direction of reaction------------------##
    reverse.index <- which(rn.info2[,4] == "reversible")
    rn.info2 <- as.data.frame(rn.info2)
    rn.info2$substrate.name[reverse.index] <- mapply(c,
      rn.info2$substrate.name[reverse.index],
      rn.info2$product.name[reverse.index])
    rn.info2$product.name[reverse.index] <- rn.info2$substrate.name[reverse.index]
    ##------------------------------------------------##
    rn.info2 <- rn.info2[,-4]
    rn <- aggregate(rn.info2, list(unlist(rn.info2[,1])),
      compose(c,unlist,unique,unlist)) %>%
      extract(-1)
    return(rn)
  }
  
  ##-------------------------------------------------------------##
  path <- sprintf("%s/%s", path, ".RevEcoR")
  if (!exists(path)) 
    dir.create(path, showWarnings = FALSE)
  metabolic.data <- sprintf("%s/%s", path, "MetabolicData.rda")
  if (!file.exists(metabolic.data)){
    message("No local metabolic data has been saved...")
    if (!refresh)
      stop("There is no local data, refresh must be TRUE, 
      ...")
    rn <- refreshData(org.pathway)
    MetabolicData <- list(rn)
    names(MetabolicData) <- org
    save(MetabolicData, file = metabolic.data)
    return(rn)
  }else{
    load(metabolic.data)
    if (refresh){
      MetabolicData[[org]]  <- refreshData(org.pathway)
      save(MetabolicData, file = metabolic.data)  
    }
    return(MetabolicData[[org]])
  }
}