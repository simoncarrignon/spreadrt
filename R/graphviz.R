#' transfrom a cascade made of numeric index into character
#' if not igraph will fill missing vertices with fake one
#' @param cascade a cascade
#' @param i an integer to identify the cascade, 0 if nothing is given
#' @return an edgelist formated to be read by `graph_from_edgelist`
formatEdgeList <- function(cascade,i=0){
    if(is.null(dim(cascade)))
        t(paste0("c",i,"-",cascade[1:2]))
    else
        t(apply(cascade,1,function(r)paste0("c",i,"-",r[1:2])))
}

