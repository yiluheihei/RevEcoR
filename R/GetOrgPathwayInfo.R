#' Get organism pathway reaction info
#' 
#' This function helps us  to obtain the specific-organism pathway map, prasing the maps to get reaction information.
#' 
#' @param org, three characters, the kegg organism code, e.g. "buc"
#' 
#' @return a three length df, consists of  reaction name, substrates and products

GetOrgPathwayInfo <- function(org){
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
  ##-------------------------------------------------------------##

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
}