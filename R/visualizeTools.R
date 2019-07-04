#' alpha
#'
#' my own implementation to add alpha to colors
alpha <- function(color,alpha){
    rgb(t(col2rgb(color)/255),alpha=alpha)
}


getCCFD <- function(){
    load("hundreds.bin")

    ca=hundreds$cascades
    res=cbind.data.frame(hundreds$summary,getScoreCascade(ca))
    png("CCDF_size.png")
    plot(result$size,(1-ecdf(result$size)(result$size))*100,log="xy",ylab="CCDF(%)",xlab="Cascade Size" )
    dev.off()
}


#'plot the CCFD following \emph{Vosoughi et al. (2018)} on a log-log axis and color it given utility
#'@param s : a two dimension vector with size of cascade and utility of the cascade for a simulation
#'@param ... : parameters passed to plot
#'@return nothing
plotCCFD_U <- function(s,pch=20,cex=2,col=NULL,legend=TRUE,...){
    s=s[order(s$size),]
    if(is.null(col)){
    br=seq(0,1,length.out=10)
    cls=cut(s$U,breaks=br)
    cols=alpha(heat.colors(length(br)),.4)
    names(cols)=levels(cls)
    col=cols[cls]
    }
    plot(c(0,s$size),c(1,(1-ecdf(s$size)(s$size)))*100,log="xy",ylab="CCDF(%)",xlab="Cascade Size",col=col,pch=pch,cex=cex,...)
    if(legend)legend("bottomleft",legend=round(br,digit=1),col=cols,pch=pch,title="U")
}

#'plot the CCFD following \emph{Vosoughi et al. (2018)} on a log-log axis and color it given utility
#'@param s : a two dimension vector with size of cascade and utility of the cascade for a simulation
#'@param ... : parameters passed to plot
#'@return nothing
pointsCCFD_U <- function(s,pch=20,cex=2,col=NULL,lenged=TRUE,...){
    s=s[order(s$size),]
    if(is.null(col)){
        br=seq(0,1,length.out=10)
        cls=cut(s$U,breaks=br)
        cols=alpha(heat.colors(length(br)),.4)
        names(cols)=levels(cls)
        col=cols[cls]
    }
    points(s$size,(1-ecdf(s$size)(s$size))*100,col=col,pch=pch,cex=cex,...)
}

#'plot the CCFD following \emph{Vosoughi et al. (2018)} on a log-log axis
#'@param s : sample of sizes
#'@param ... : parameters passed to plot
#'@return nothing
plotCCFD <- function(s,pch=20,cex=1,xlab="",ylab="CCDF(%)",...){
    ccfd=ccfd(s)
    plot(ccfd[,"x"],ccfd[,"y"],log="xy",,pch=pch,ylab=ylab,cex=cex,xlab=xlab,...)
}

#' add points of the CCFD following \emph{Vosoughi et al. (2018)}
#'@param s : sample of sizes
#'@param ... : parameters passed to plot
#'@return nothing
pointsCCFD <- function(s,pch=20,cex=1,...){
    ccfd=ccfd(s)
    points(ccfd[,"x"],ccfd[,"y"],pch=pch,cex=cex,...)
}



#' Calculate the CCFD following \emph{Vosoughi et al. (2018)}
#'@param size : sample of sizes
#'@return a order table with two column with first column : frequency and second column the probality of the frequency  
ccfd <- function(size){
    total=length(size) 
    size=size[order(size)]
    counts=unique(size)
    id=unique(match(counts,size))
    p=sapply(id,function(i)length(size[i:total])) 
    p=p/total * 100
    x=counts
    y=p
    return(cbind(x,y))
}


#' return a vector of color that follow the values of utilities 
#'@param utilites : a vector with value we want to match
#'@return a vector of length = `length(utilities)` of colors with names=sort(names(utilities))
getUCols <- function(utilities){
    unique_u = unique(utilities) #all possible utility as defined in the default function `cascades3D`
    cols=heat.colors(length(unique_u))
    names(cols)=sort(unique_u)
    return(cols[as.character(utilities)])
}



#' plot posteriors distribution against priors
#'@param A : a vector with posterior
#'@param B : a vector with posterior 
#'@param prior : a vector with posterior 
plot2dens <- function(A=NULL,B=NULL,C=NULL,from=NULL,to=NULL,prior=NULL,cols=c(alpha("red",.8),alpha("blue",.8),alpha("yellow",.8)),hdr=F,yaxt=NULL,...){

    denseP=NULL
    denseA=NULL
    denseB=NULL
    denseC=NULL
    if(!is.null(prior))prior=prior[!is.na(prior)]
    if(is.null(yaxt))yaxt="n"
    if(is.null(from))from=min(A,B,prior)
    if(is.null(to))to=max(A,B,prior)
    if(!is.null(A))denseA=density(A,from=from,to=to)
    if(!is.null(B))denseB=density(B,from=from,to=to)
    if(!is.null(C))denseC=density(C,from=from,to=to)
    if(length(prior)==2)denseP=density(runif(100000,prior[1],prior[2]),from=from,to=to)
    else if(!is.null(prior))denseP=density(prior,from=from,to=to)
    print("donde")

    if(is.null(names(cols)))names(cols)=c("P","A","B")
    rangex=range(denseB$x,denseA$x,denseP$x,denseC$x)
    rangey=range(0,denseB$y,denseA$y,denseP$y,denseC$y)
    stepy=max(rangey)*0.2
    if(hdr)miny=-1*stepy else miny=0
    plot(denseA,ylim=c(miny,max(rangey)),xlim=rangex,type="n",xaxt="n",yaxt=yaxt,ylab="",...)
    axis(1)
    mtext("Density",2,1)
    if(!is.null(prior))
        polygon(c(from,denseP$x,to),c(0,denseP$y,0),col=cols["P"],lwd=2)
    if(!is.null(C))
        polygon(c(from,denseC$x,to),c(0,denseC$y,0),col=cols["C"],lwd=2)
    if(!is.null(A)){
        #polygon(c(from,denseA$x,to),c(0,denseA$y,0),col="white",border=NA,lwd=2)
        if(hdr){
            hdstaA=hdr(A,prob=c(75,95),lambda=0.9)
            hdrA=hdstaA$hdr
            print(hdrA)
            polygon(c(hdrA[1,1],hdrA[1,1],hdrA[1,2],hdrA[1,2]),c(-.2*stepy,-.9*stepy,-.9*stepy,-.2*stepy),col=cols["A"],lwd=1)

            polygon(c(hdrA[2,1],hdrA[2,1],hdrA[2,2],hdrA[2,2]),c(-.2*stepy,-.9*stepy,-.9*stepy,-.2*stepy),col=cols["A"],lwd=1)
            polygon(c(hdrA[2,1],hdrA[2,1],hdrA[2,2],hdrA[2,2]),c(-.2*stepy,-.9*stepy,-.9*stepy,-.2*stepy),col=cols["A"],lwd=1)
            segments(hdstaA$mode,-.3*stepy,hdstaA$mode,-.8*stepy)
            middle=(-.2*stepy+-.9*stepy)/2
            segments(min(A),middle,hdrA[1,1],middle)
            segments(hdrA[1,2],middle,max(A),middle)
            polygon(c(from,denseA$x,to),c(0,denseA$y,0),col=cols["A"],lwd=2)
        }
        else
            polygon(c(from,denseA$x,to),c(0,denseA$y,0),col=cols["A"],lwd=2)
    }
    if(!is.null(B))
        polygon(c(from,denseB$x,to),c(0,denseB$y,0),col=cols["B"],lwd=2)


}
