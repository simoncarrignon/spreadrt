
#' getMaxBreadth: 
#' Gives the max breadth of a tree
#' @param cascade: a edgelist that represent a cascade
#' @return integer 
getMaxBreadth <- function(cascade){
    if(nrow(cascade)<3) return(1)
    root=cascade[1,1]
    cascade=cascade[-1,]
    bmax=1
    while(!is.null(root)){

        child=cascade[cascade[,1] %in% root,] #get all child of current nodes

        if(is.null(dim(child)) & length(child==3))  #if there is just one child child will be a dimension less vector of size 3, with child[2] the id of the child
            root = child[2]
        else if(!is.null(dim(child))){
            if(dim(child)[1] == 0) #if there is no mor child we will have an empty table with three columns: we reach all the leaves
                return(bmax)
            else{
                root=child[,2]
                b=nrow(child)
                if(b>bmax)bmax=b
            }
        }

    }

}
#' getDepth: 
#' Gives the depth of a tree
#' @param cascade: a edgelist that represent a cascade
#' @return integer 
#' @export getDepth 
getDepth <- function(cascade){
    if(nrow(cascade)==1)
        return(1)
    if(is.null(dim(cascade[-1,])))
        return(2)
    max(unlist(depths(cascade[-1,],cascade[1,1])))
}

getDepth <- function(cascade){
    if(nrow(cascade)==1)
        return(1)
    if(is.null(dim(cascade[-1,])))
        return(2)
    max(unlist(depths(cascade[-1,],cascade[1,1])))
}


#' depths: 
#' the recursive function called to compute the depth
#' @param cascade: a edgelist that represent a cascade
#' @param root: the id of the root node
#' @return list of depth of each branch
#' @export depths 
depths <- function(cascade,root){
    childs=cascade[cascade[,1] == root,2]
    if(length(childs)==0)
        return(1) 
    else {
        sapply(childs,function(c)return(1+max(depths(cascade,c))))
    }
}

#' Binned Frequencies distance
#' distance between the frequencies of a and b, the frequencies beings calculates for logarithmics categories
#'@param a and b: sample of any size 
#'@return a list of dist
#' @export euclidediff 
euclidediff <- function(a,b,nbin=20) {
    m=max(a,b)
    pm=round(log10(m))+1
    bins=10^(seq(0,pm,length.out = nbin))
    taba=prop.table(table(cut(a,breaks=bins,include.lowest = T)))
    tabb=prop.table(table(cut(b,breaks=bins,include.lowest = T)))
    sum((taba-tabb)^2)
}

#' Stupid Diff Function
#' distance between two sample a and b
#'@param a and b: sample of any size 
#'@return a list of dist
diff <- function(a,b) {
    smax=max(a,b)
    res=abs(log(ecdf(a)(1:smax))-log(ecdf(b)(1:smax)))
    res[!is.infinite(res)]

}

#' Quantil Diff Function
#' distance between two samples a and b by looking at the difference between the log o all percentile of two sample
#'@param a and b: sample of any size 
#'@return a number
quantilediff <- function(a,b,log=T,probs=seq(0,1,.01)) {
    sqrt(sum((quantile(log(b),probs=probs)-quantile(log(a),probs=probs))^2))
}

#' kindof KL
#' distance between two samples a and b by looking at the difference between the log o all percentile of two sample
#'@param a and b: sample of any size 
#'@return a number
kokldiff <- function(a,b,probs=seq(0,1,.01)) {
    sum(log(quantile(a,probs=probs)/quantile(b,probs=probs))*quantile(a,probs=probs))
}


#' More kindof KL
#' distance between two samples a and b by looking at the difference between the log o all percentile of two sample
#'@param a and b: sample of any size 
#'@return a number
kldiff <- function(a,b) {
    rs=range(a,b) #this is to find the range where our data are defined, will be used to avoid

    x=rs[1]:rs[2]  #generate the space where both function are defined

    ### Kind of Laplace Smoothing 
    a=c(a,x)
    b=c(b,x)

    #to calculate the KL-Divergence we need to compare the probability function and not the sample
    P=approxfun(density(a,from=rs[1],to=rs[2]))  #we need then first to estimate thoses probabilty function of both sample (the data and the simulation). To do so we first get the probability density function using a kernel regression and from this interpolate a linear function
    Q=approxfun(density(b,from=rs[1],to=rs[2]))
    Px=P(x) 
    Qx=Q(x)
    ndef=c(which(Qx==0),which(Px==0)) #with the laplace smoothing that should'n make any dif
    if(length(ndef)>0)x=x[-ndef]  #we exclude the part of X where P(x) or Q(x) is not defined. This may not be right but it's a good wait to avoid the sum to be infinite or not defined.
    sum(log(P(x)/Q(x))*P(x))
}


tweetvsUnB <- function(mat,u,betamean){
    cols=topo.colors(length(betamean))
    names(cols)=betamean
    sapply(as.character(betamean),function(m)lines(u,mat[,m],col=alpha(cols[as.character(m)],.8),lwd=3))
}

#' probaSelection a return the probability for a list of tweet of given  utility to be retweeted given a value of beta
#' @param u : vector of utility of tweet
#' @param beta : value of beta 
#' @return a vector of probability of the size length(u) 
probaSelection <- function(u,beta){
            p_i=exp(u * beta)
            p_i=p_i/sum(p_i)
            return(p_i)
}

#' generateCascades initialize a list of cascades of retweet with only the first tweet
#' @param nc : number of cascades
#' @param seeders : list or number of possible agents that will start the cascade
#' @param random : list or number of possible rumors where to choose the cascades
#' @param time : the time at wich start the cascades
#' @return a list of size nc with the id of the first agents and the rumor
generateCascades <- function(nc,seeders,rumours,time){
    stopifnot(seeders*rumours>=nc) #check if we can theoretically generate the good number of initial cascade

    #there are probably better way to do than (by choosing unique pair of rumor/seeder insteed of randomly choose and check afterward for duplicate)
    nseeders=sample.int(seeders,nc,replace=T)
    nrumors=sample.int(rumours,nc,replace=T)
    ncascades_id=data.frame(seeder=nseeders,rumor=nrumors)
    ncascades_id=unique(ncascades_id)
    while(nrow(ncascades_id)<nc){
        nseeders=sample.int(seeders,1)
        nrumors=sample.int(rumours,1)
        ncascades_id=rbind.data.frame(ncascades_id,c(seeder=nseeders,rumor=nrumors))
        ncascades_id=unique(ncascades_id)
    }
    ncascades=lapply(1:nc,function(i)array(c(ncascades_id$seeder[i],ncascades_id$seeder[i],time),c(1,3)))
    return(list(cascades=ncascades,cascades_id=ncascades_id))
}

#' tomat from a data frame create a matrix
#' @param dt : dataframe
#' @param fun : a function returning a single number from a sample
#' @param measure : the measurement on which will be applied the function fun
#' @return a matrix
tomat <- function(dt,fun=mean,measure="size")tapply(dt[,measure],dt[,c("U","beta")],fun)

#' removeDuplicate return index of the vector v without the replicate and by randomly selecting amoung the replicate
#' @param v : a vector
#' @return a vector of indice of v
removeDuplicate <- function(v){
    id=1:length(v)
    nid=c()
    for(i in unique(v)){
        cnt=length(v[v==i])
        if(cnt>1)nid=c(nid,sample.int(id[v==i],1))
        else nid=c(nid,id[v==i])
    }
    return(nid)
}


#' getSummaryCascade: from a list of cascade, return a list with the size, depth and max breadth of each cascades
#' @param cascade : a list of cascades
#' @return a  table with the column define in metrics, ("size","depth","breadth") by default, and a number of row equal to the length of the input list
getSummaryCascade <- function(cascade,metrics=c("size","depth","breadth")){
    listr=sapply(1:length(cascade),function(c)
                 {
                     one=cascade[[c]]
                     res=c()
                     if(sum(metrics=="size")>=1)res["size"]=nrow(one)
                     if(sum(metrics=="breadth")>=1)res["breadth"]=getMaxBreadth(one)
                     if(sum(metrics=="depth")>=1)res["depth"]=getDepth(one)
                     return(res)
                 }
    )
    if(is.null(dim(listr))){
        dim(listr)=c(length(cascade),1)
        colnames(listr)=metrics
        return(listr)
    }
    return(t(listr))
}


#' Generate a random three partition 
threePartition <- function(n){
    a=runif(n,0,1)
    b=runif(n,0,1-a)
    c=1-(a+b)
    res=cbind(a,b,c)
   return(t(apply(res,1,function(r)r[sample(3)])))
}

#' Generate a table  n subvector of random size 
#'@param n number of row
#'@param d number of column
#'@return a table of dim nxd
generalPartition <- function(n,d){
    res=array(dim=c(n,0))
    #a=runif(n,0,1)
    #res=cbind(res,a)
    #a=runif(n,0,1-res)
    #res=cbind(res,a)
    for(i in 1:(d-1)){
        new=runif(n,0,1-apply(res,1,sum))
        res=cbind(res,new)
    }
    new=1-apply(res,1,sum)
    res=cbind(res,new)
    return(t(apply(res,1,function(r)r[sample(d)])))
}


#' Generate a Distribution of n value splited in   reparted among all class of the vector `cls` following the frequenceies `freq`
#'@param cls all possible classe
#'@param n total number of value
#'@param freq final frequencies for each class
#'@return a vector of size `n` with value from `cls` that respect the frequencies in `freq`
generateDistribeFromFrequencies <- function(cls,n,freq){
    u=rep(cls,freq*n)
    tmp_n=sum(table(u))
    if(tmp_n<n){ ##In cas of missing value we ad to ad more
        miss=n-tmp_n
        u=c(u,sample(cls,miss,repl=T))
        u
    }
    else u
}


plotDistributGivenScore <- function(listoffullscores,listofscore,m,rumors=F,...){
    if(rumors && m !="size")return
    names(listoffullscores)=listofscore
    cols=alpha(topo.colors(length(listofscore)),.5)
    names(cols)=as.character(sort(listofscore))
   if(rumors) plotCCFD(allru$size,cex=1,pch=1,...)
   else plotCCFD(allca[[m]],cex=1,pch=1,...)
    n=lapply(1:length(listoffullscores),function(l)
             if(rumors)pointsCCFD(getRumors(listoffullscores[[l]]$rd$rd),col=cols[names(listoffullscores)[[l]]],type="l",lwd=3)
             else pointsCCFD(listoffullscores[[l]]$rd$rd[[m]],col=cols[names(listoffullscores)[[l]]],type="l",lwd=3)
             )
   if(rumors) pointsCCFD(allru$size,cex=1,pch=1)
   else pointsCCFD(allca[[m]],cex=1,pch=1)
}

#' add the size of the cascade of the same rumors togethers
getRumors <- function(tablres,metric="size")tapply(tablres[[metric]],tablres$rumor,sum)

#' simple function to handle parameter space made of multiple dimension parameter
#'@param vec a vector potentially multidimensional
#'@param filter the index of the vector that we want to keep
filter <- function(vec,filter)if(is.null(dim(vec)))vec[filter] else if(length(dim(vec))==2)vec[filter,]
