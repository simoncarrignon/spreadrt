library(devtools)
load_all("../spread")


getAllScoresFromFolder <- function(){
    expefold="correctmu/"
    #load(paste0(file.path(expefold,basename(expefold)),"/parameters.bin"))
    totest=c(1:94)
    allrumorscores=sapply(totest,function(i){
                          load(file.path(expefold,paste0(basename(expefold),i),"scores.bin"))
                          sapply(scores,function(l)l[["rumors"]]) 
})

    allcascadescores=sapply(totest,function(i){
                            load(file.path(expefold,paste0(basename(expefold),i),"scores.bin"))
                            sapply(scores,function(l)l[["cascades"]]) 
})

    allparameters=lapply(totest,function(i){
                         load(file.path(expefold,paste0(basename(expefold),i),"parameters.bin"))
                         as.data.frame(parameters)
})
    allparameters.dataframe=do.call("rbind.data.frame",allparameters) 
}

#Those are the important and full results
getTheAllStuffPreviouslyAlreadyExtracted <- function(){
    load("allparameters.df") 
    load("allscores.cascade") 
    load("allscores.rumors") 
}


generateResultsGraphRandom <- function(){
    png("density_bothscore.png",width=2*480)
    par(mfrow=c(1,2),mar=c(4,4,2,1))
    plot(density(allcascadescores),main="size of rumors",xlab="distance to data")
    plot(density(allrumorscores),main="size of cascades",xlab="distance to data")
    dev.off()

    png("correlation_bothscore.png")
    par(mfrow=c(1,1),mar=c(4,4,2,1))
    plot(allcascadescores ~ allrumorscores,col=alpha("black",.2),pch=20)
    dev.off()

    ###echec:
    rscors=sapply(backup_cascades,function(l)l$scores$rumors)
    rscors.c=sapply(backup_cascades,function(l)l$scores$cascades)
    cols=alpha(topo.colors(length(rscors)),.5)
    names(cols)=as.character(sort(rscors))
    plotCCFD(allru$size,cex=1,pch=1,main="simu vs rumor size")
    n=lapply(backup_cascades,function(l)
             pointsCCFD(l$rd$size,col=cols[as.character(l$scores$rumors)],type="l",lwd=3)
             )
    ####fin de lechec



    # use the same syntaxe that in the vignet:
    posteriors.r=allparameters.dataframe[rank(allrumorscores)<500,]
    posteriors.c=allparameters.dataframe[rank(allcascadescores)<500,]


    png("rumors_postcheck.png")
    plotCCFD(allru$size)
    nan=apply(posteriors.r,1,function(u){pointsCCFD(randomCascades(Nmin=u["Nmin"],Nmax=u["Nmax"],mu=u["mu"],t_step=u["t_step"],tau=u["tau"])$size,type="l",col=alpha("chartreuse",.4))})
    dev.off()

    png("rumors_median_postcheck.png")
    median_best_param=apply(posteriors.r,2,median)
    plotCCFD(allru$size)
    na=replicate(100,pointsCCFD(randomCascades(Nmin=median_best_param["Nmin"],Nmax=median_best_param["Nmax"],mu=median_best_param["mu"],t_step=median_best_param["t_step"],tau=median_best_param["tau"])$size,type="l",col=alpha("red",.1)))
    dev.off()

    png("rumors_best_postcheck.png")
    best_param=unlist(allparameters.dataframe[which(rank(allrumorscores) == 1),])
    plotCCFD(allru$size)
    na=replicate(100,pointsCCFD(randomCascades(Nmin=best_param["Nmin"],Nmax=best_param["Nmax"],mu=best_param["mu"],t_step=best_param["t_step"],tau=best_param["tau"])$size,type="l",col=alpha("red",.1)))
    dev.off()


    png("cascades_postcheck.png")
    plotCCFD(allca$size)
    nan=apply(posteriors.c,1,function(u){pointsCCFD(randomCascades(Nmin=u["Nmin"],Nmax=u["Nmax"],mu=u["mu"],t_step=u["t_step"],tau=u["tau"])$size,type="l",col=alpha("chartreuse",.4))})
    dev.off()


    png("cascades_median_postcheck.png")
    median_best_param=apply(posteriors.c,2,median)
    plotCCFD(allca$size)
    na=replicate(100,pointsCCFD(randomCascades(Nmin=median_best_param["Nmin"],Nmax=median_best_param["Nmax"],mu=median_best_param["mu"],t_step=median_best_param["t_step"],tau=median_best_param["tau"])$size,type="l",col=alpha("red",.1)))
    dev.off()

    png("cascades_best_postcheck.png")
    best_param=unlist(allparameters.dataframe[which(rank(allcascadescores) == 1),])
    plotCCFD(allca$size)
    na=replicate(100,pointsCCFD(randomCascades(Nmin=best_param["Nmin"],Nmax=best_param["Nmax"],mu=best_param["mu"],t_step=best_param["t_step"],tau=best_param["tau"])$size,type="l",col=alpha("red",.1)))
    dev.off()

    png("posterior_distrib.png",width=450,height=3*330)
    par(mfrow=c(3,2))
    sapply(colnames(allparameters.dataframe),function(p)plot2dens(A=allparameters.dataframe[rank(allrumorscores)<500,p],B=allparameters.dataframe[rank(allcascadescores)<500,p],prior = allparameters.dataframe[,p],main=p,xlab=""))
    plot(1,1,type="n",axes=F)
    legend("center",legend=c("prior","cascades","rumors"),fill=c("red","yellow","blue"))
    dev.off()

    png("mu_distrib.png",width=480,height=320)
    par(mfrow=c(1,2))
    plot(density(allparameters.dataframe$mu[rank(allrumorscores)<500],from=0.00001),main="rumors",xlab=expression(mu))
    plot(density(allparameters.dataframe$mu[rank(allcascadescores)<500],from=0.00001),main="cascades",xlab=expression(mu))
    dev.off()

    png("param_corel_rumor.png",width=480)
    par(mfrow=c(1,2))
    plot(allparameters.dataframe[rank(allrumorscores)<500,])
    dev.off()
    png("param_corel_ca.png",width=480)
    plot(allparameters.dataframe[rank(allcascadescores)<500,])
    dev.off()

}

generateResultsGraph3D <- function(){
    library(devtools)
    load_all("../spread/")
    #in this case the data are in : sshfs dlmn:/gpfs/scratch/bsc21/bsc21394/twitter-spread/abcdir/laplace laplace/
    mainfoldresults="../laplace/"
    #With all the folder we can now get all the scores
    allscores = getAllscores(mainfoldresults,lim=0:10)
    allparameters.dataframe = getAllparameters(mainfoldresults,lim=0:10)

    #With all the folder we can now get all the scores
    foldrum="../testneutralRumors"
    rumors.scores = getAllscores(foldrum,idscores = c("qc"),metrics=c("false","true","size"),lim=0:6)
    rumors.parameters = getAllparameters(foldrum,lim=0:6)
    rumors.posteriors=getAllposteriors(rumors.scores,rumors.parameters,500)
    rumors.posteriors.scores.kc=rumors.scores$kc$size[rank(rumors.scores$kc$size,ties="first")<500]
    rumors.posteriors.scores.qc=rumors.scores$qc$size[rank(rumors.scores$qc$size,ties="first")<500]

    foldrum="../testConf"
    neutmod.scores = getAllscores(foldrum,idscores = c("kc","qc"),metrics="size",lim=2:101)
    neutmod.parameters = getAllparameters(foldrum,log=T,lim=2:101)
    neutmod.posteriors=getAllposteriors(neutmod.scores,neutmod.parameters,500)
    neutmod.posteriors.scores.kc=neutmod.scores$kc$size[rank(neutmod.scores$kc$size,ties="first")<500]
    neutmod.posteriors.scores.qc=neutmod.scores$qc$size[rank(neutmod.scores$qc$size,ties="first")<500]


    #once we have all the scrore and all the parameter we can start playing with
    png("density_scores_3D.png",height=3*480,width=3*480) 
    par(mfrow=c(5,3),mar=c(4,4,2,1))
    lapply(idscores,function(s)
           lapply(metrics,function(m)
                  {
                      cur=allscores[[s]][m]
                      cur=cur[!is.na(cur)]
                      try(plot(density(cur),main=paste(m," of cascades"),xlab=s))
                  }
           )
           )
    dev.off()

    png("correlation_scores_3D.png",height=3*480,width=8*480)
    par(mfrow=c(3,10),mar=c(4,4,2,1))
           lapply(metrics,function(m){
                  try(plot( allscores[["kc"]][[m]] ~ allscores[["qc"]][[m]]   ,log="xy",main=m,xlab="qc",ylab="kc",col=alpha("black",.2),pch=20))
                  try(plot( allscores[["kc"]][[m]] ~ allscores[["kokc"]][[m]] ,log="xy",main=m,xlab="kokc",ylab="kc",col=alpha("black",.2),pch=20))
                  try(plot( allscores[["qc"]][[m]] ~ allscores[["kokc"]][[m]] ,log="xy",main=m,xlab="kokc",ylab="qc",col=alpha("black",.2),pch=20))
                  try(plot( allscores[["kc"]][[m]] ~ allscores[["kokcrev"]][[m]] ,log="xy",main=m,xlab="kokcrev",ylab="kc",col=alpha("black",.2),pch=20))
                  try(plot( allscores[["qc"]][[m]] ~ allscores[["kokcrev"]][[m]] ,log="xy",main=m,xlab="kokcrev",ylab="qc",col=alpha("black",.2),pch=20))
                  try(plot( allscores[["kokc"]][[m]] ~ allscores[["kokcrev"]][[m]] ,log="xy",main=m,xlab="kokcrev",ylab="kokc",col=alpha("black",.2),pch=20))
                  try(plot( allscores[["kcrev"]][[m]] ~ allscores[["qc"]][[m]]   ,log="xy",main=m,xlab="qc",ylab="kcrev",col=alpha("black",.2),pch=20))
                  try(plot( allscores[["kcrev"]][[m]] ~ allscores[["kokc"]][[m]] ,log="xy",main=m,xlab="kokc",ylab="kcrev",col=alpha("black",.2),pch=20))
                  try(plot( allscores[["kcrev"]][[m]] ~ allscores[["kokcrev"]][[m]] ,log="xy",main=m,xlab="kokcrev",ylab="kcrev",col=alpha("black",.2),pch=20))
                  try(plot( allscores[["kcrev"]][[m]] ~ allscores[["kc"]][[m]] ,log="xy",main=m,xlab="kc",ylab="kcrev",col=alpha("black",.2),pch=20))
           })
    dev.off()


    posteriors=lapply(allscores,lapply,function(u)allparameters.dataframe[rank(u,ties="first")<500,])
    posteriorsb=getAllposteriors(allscores,allparameters.dataframe,500)

    getbest(allscores = rumors.scores,allparameters.dataframe = rumors.parameters,"qc","size")
    #as calculate all that is heavy and costly, better keeping it somewhere safe 
    #save(file="LAPLACEposteriors3D5Scores3Measurments.bin",posteriors)
    #save(file="LAPLACElastallscores3D5Scores3Measurments.bin",allscores)
    #save(file="LAPLACEparam.lastallscores3D5Scores3Measurments.bin",allparameters.dataframe)


    for(s in idscores){
        for(m in metrics){
            png(width=3*480,height=3*480,paste0("parametercorrelation_",s,"_",m,".png"))
            plot(posteriors[[s]][[m]])
            dev.off()
        }
    }
    for(param in colnames(allparameters.dataframe)){
        #dev.new()
        png(paste0("all_posteriors_",param,".png"),width=4*480,height=4*480)
        par(mfrow=c(5,3),mar=c(5,5,1,1))
        for(s in idscores){
            posteriors.s=posteriors[[s]]
            for(m in metrics){
                posteriors.qc = posteriors.s[[m]]
                prior=density(allparameters.dataframe[,param]) #prior distribution
                postC=density(posteriors.qc[,param])  #posterior distribution for distance to cascade size
                rangex=range(prior$x,postC$x)
                rangey=range(prior$y,postC$y)
                plot(prior,ylim=rangey,xlim=rangex,lwd=5,main=paste("scores=",s,"m=",m,"param",param),xlab=param)
                lines(postC,col="red",lwd=5)
            }
        }
        dev.off()
    }
    sapply(dev.list(),dev.off)

    ### POSTERIORS CHECKS
util=grep("utility.*",colnames(allparameters.dataframe))
betad=grep("betaDistrib.*",colnames(allparameters.dataframe))

    png("best_cascade_postcheck_3D.png",width=2*480,height=3*480)
    par(mfrow=c(5,3))
    lapply(names(allscores),function(s)lapply(names(allscores[[s]]),function(m){
                                              print(paste(s,m))
                                              best_param=unlist(allparameters.dataframe[which(rank(allscores[[s]][[m]],ties="first") == 1),])
                                              plotCCFD(allca[[m]],main=paste(m,s))
                                              na=replicate(50,pointsCCFD(
                                                                          cascades3D(
                                                                                     log=F,
                                                                                     N=best_param["N"],
                                                                                     stime=best_param["stime"],
                                                                                     dtime=best_param["dtime"],
                                                                                     IC=best_param["IC"],
                                                                                     Nmax=best_param["Nmax"],
                                                                                     R=best_param["R"],
                                                                                     time=best_param["repetition"],
                                                                                     captl=best_param["captl"],
                                                                                     utility = generateDistribeFromFrequencies(prior_utility,best_param["R"],best_param[util]),
                                                                                     betadistrib = generateDistribeFromFrequencies(prior_beta,best_param["N"],best_param[betad])
                                                                                     )[[m]],type="l",col=alpha("red",.1)))
                  }))
    dev.off()

    onbest2time=cascades3D(
                      log=F,
                      N=best_param["N"]*100,
                      stime=best_param["stime"],
                      dtime=best_param["dtime"],
                      IC=best_param["IC"],
                      Nmax=best_param["Nmax"],
                      R=best_param["R"],
                      time=best_param["repetition"]*2,
                      captl=best_param["captl"],
                      utility = generateDistribeFromFrequencies(prior_utility,best_param["R"],best_param[util]),
                      betadistrib = generateDistribeFromFrequencies(prior_beta,best_param["N"]*100,best_param[betad])
                      )
    onerand=cascades3D(
                      log=F,
                      N=best_param["N"],
                      stime=best_param["stime"],
                      dtime=best_param["dtime"],
                      IC=best_param["IC"],
                      Nmax=best_param["Nmax"],
                      R=best_param["R"],
                      time=best_param["repetition"],
                      captl=best_param["captl"],
                      utility = generateDistribeFromFrequencies(c(-1,0,1),best_param["R"],best_param[util]),
                      betadistrib = generateDistribeFromFrequencies(c(0,0,0,0,0),best_param["N"],best_param[betad])
                      )

    par(mfrow=c(1,2))
    plotCCFD(allca$size,main="Cascades size",axes=F)
    pointsCCFD(mixedca$size,col="orange")
    pointsCCFD(falseca$size,col="red")  
    pointsCCFD(trueca$size ,col="green")
    #pointsCCFD(onbest$size,pch=4)
    pointsCCFD(onbest$size[onbest$U==0],col="orange",pch=4)
    pointsCCFD(onbest$size[onbest$U<0] ,col="red",pch=4)
    pointsCCFD(onbest$size[onbest$U>0] ,col="green",pch=4)
 legend("bottomleft",pch=c(20,20,20,1,4),col=c("red","green","orange","black","black"),legend=c("false (U = -1)","true(U = 1)","mixed (U = 0)","data","model") )
    plotCCFD(allru$size,type="n",main="Rumor size")
    pointsCCFD(mixedru$size,col="orange")
    pointsCCFD(falseru$size,col="red")  
    pointsCCFD(trueru$size ,col="green")
    #pointsCCFD(onbest$size,pch=4)
    pointsCCFD(onbest2time$size[onbest$U==0],col="orange",pch=4)
    pointsCCFD(onbest2time$size[onbest$U<0] ,col="red",pch=4)
    pointsCCFD(onbest2time$size[onbest$U>0] ,col="green",pch=4)
 legend("bottomleft",pch=c(20,20,20,1,4),col=c("red","green","orange","black","black"),legend=c("false (U = -1)","true(U = 1)","mixed (U = 0)","data","model") )

    plotCCFD(allca$size,type="n")
    pointsCCFD(mixedca$size,col="orange")
    pointsCCFD(falseca$size,col="red")  
    pointsCCFD(trueca$size ,col="green")
    pointsCCFD(onerand$size[onerand$U==0],col="orange",pch=4)
    pointsCCFD(onerand$size[onerand$U<0],col="red",pch=4)
    pointsCCFD(onerand$size[onerand$U>0],col="green",pch=4)

    png("all_cascade_postcheck_3D.png",height=2*480,width=1.5*480)
    par(mfrow=c(5,3))
    lapply(names(allscores),function(s)lapply(names(allscores[[s]]),function(m){
                                              print(paste("all poster",s,m))
                                              plotCCFD(allca[[m]])
                                              apply(posteriors[[s]][[m]],1,function(best_param)
                                                    pointsCCFD(
                                                               cascades3D(
                                                                          log=F,
                                                                          N=best_param["N"],
                                                                          stime=best_param["stime"],
                                                                          dtime=best_param["dtime"],
                                                                          IC=best_param["IC"],
                                                                          Nmax=best_param["Nmax"],
                                                                          R=best_param["R"],
                                                                          time=best_param["repetition"],
                                                                          captl=best_param["captl"],
                                                                          utility = generateDistribeFromFrequencies(prior_utility,best_param["R"],best_param[util]),
                                                                          betadistrib = generateDistribeFromFrequencies(prior_beta,best_param["N"],best_param[betad])
                                                                          )[[m]],type="l",col=alpha("red",.1)))
                  }))
    dev.off()


    png("median_cascade_postcheck_3D.png",width=2*480,height=3*480)
    par(mfrow=c(5,3))
    lapply(names(allscores),function(s)lapply(names(allscores[[s]]),function(m){
                                              print(paste("median",s,m))
                                              median_best_param=apply(posteriors[[s]][[m]],2,median)
                                              plotCCFD(allca[[m]],main=paste(m,s))
                                              na=replicate(50,pointsCCFD(
                                                                          cascades3D(
                                                                                     log=F,
                                                                                     N=median_best_param["N"],
                                                                                     stime=median_best_param["stime"],
                                                                                     dtime=median_best_param["dtime"],
                                                                                     IC=median_best_param["IC"],
                                                                                     Nmax=median_best_param["Nmax"],
                                                                                     R=median_best_param["R"],
                                                                                     time=median_best_param["repetition"],
                                                                                     captl=median_best_param["captl"],
                                                                                     utility = generateDistribeFromFrequencies(prior_utility,median_best_param["R"],median_best_param[util]),
                                                                                     betadistrib = generateDistribeFromFrequencies(prior_beta,median_best_param["N"],median_best_param[betad])
                                                                                     )[[m]],type="l",col=alpha("red",.1)))
                  }))
    dev.off()
    ##Calcualte the mode from S-O
    Mode <- function(x) {
          ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    png("Mode_cascade_postcheck_3D.png",width=2*480,height=3*480)
    par(mfrow=c(5,3))
    lapply(names(allscores),function(s)lapply(names(allscores[[s]]),function(m){
                                              print(paste("Mode",s,m))
                                              Mode_best_param=apply(posteriors[[s]][[m]],2,Mode)
                                              plotCCFD(allca[[m]],main=paste(m,s))
                                              na=replicate(50,pointsCCFD(
                                                                          cascades3D(
                                                                                     log=F,
                                                                                     N=Mode_best_param["N"],
                                                                                     stime=Mode_best_param["stime"],
                                                                                     dtime=Mode_best_param["dtime"],
                                                                                     IC=Mode_best_param["IC"],
                                                                                     Nmax=Mode_best_param["Nmax"],
                                                                                     R=Mode_best_param["R"],
                                                                                     time=Mode_best_param["repetition"],
                                                                                     captl=Mode_best_param["captl"],
                                                                                     utility = generateDistribeFromFrequencies(prior_utility,Mode_best_param["R"],Mode_best_param[util]),
                                                                                     betadistrib = generateDistribeFromFrequencies(prior_beta,Mode_best_param["N"],Mode_best_param[betad])
                                                                                     )[[m]],type="l",col=alpha("red",.1)))
                  }))
    dev.off()

    repbest=replicate(3,cascades3D(
                                   log=F,
                                   N=best_param["N"],
                                   stime=best_param["stime"],
                                   dtime=best_param["dtime"],
                                   IC=best_param["IC"],
                                   Nmax=best_param["Nmax"],
                                   R=best_param["R"],
                                   time=best_param["repetition"],
                                   captl=best_param["captl"],
                                   utility = generateDistribeFromFrequencies(prior_utility,best_param["R"],best_param[util]),
                                   betadistrib = generateDistribeFromFrequencies(prior_beta,best_param["N"],best_param[betad])
                                   )[["size"]])

    prm=allparameters.dataframe[1,]
    testna=replicate(100,cascades3D(
                                   log=F,
                                   N=prm["N"],
                                   stime=prm["stime"],
                                   dtime=prm["dtime"],
                                   IC=prm["IC"],
                                   Nmax=prm["Nmax"],
                                   R=prm["R"],
                                   time=prm["repetition"],
                                   captl=prm["captl"],
                                   utility = generateDistribeFromFrequencies(prior_utility,prm["R"],prm[util]),
                                   betadistrib = generateDistribeFromFrequencies(prior_beta,prm["N"],prm[betad])
                                   )
    )

    tscore=sapply(mes,function(ma)apply(testna,2,function(cl,m=ma)spread::kldiff(cl[[m]],allca[,m]))) 
    revtscore=sapply(mes,function(ma)apply(testna,2,function(cl,m=ma)spread::kldiff(allca[,m],cl[[m]]))) 


    atscore=sapply(mes,function(ma)apply(testna,2,function(cl,m=ma)spread::kokldiff(cl[[m]],allca[,m]))) 
    arevtscore=sapply(mes,function(ma)apply(testna,2,function(cl,m=ma)spread::kokldiff(allca[,m],cl[[m]]))) 

    couls=heat.colors(length(unique(tscore)))
    names(couls)=sort(unique(tscore))
    plotCCFD(allca$size)
    for(i in 1:ncol(testna)) pointsCCFD(testna[,i]$size,col=couls[as.character(tscore[i,"size"])],type="l",lwd=.5)

    png("cascades_postcheck_3D.png")
    plotCCFD(allca$size)
    nan=apply(posteriors.c,1,function(u){pointsCCFD(randomCascades(Nmin=u["Nmin"],Nmax=u["Nmax"],mu=u["mu"],t_step=u["t_step"],tau=u["tau"])$size,type="l",col=alpha("chartreuse",.4))})
    dev.off()



    png("posterior_distrib_3D.png",width=450,height=3*330)
    par(mfrow=c(3,2))
    sapply(colnames(allparameters.dataframe),function(p)plot2dens(A=allparameters.dataframe[rank(allrumorscores)<500,p],B=allparameters.dataframe[rank(allcascadescores)<500,p],prior = allparameters.dataframe[,p],main=p,xlab=""))
    plot(1,1,type="n",axes=F)
    legend("center",legend=c("prior","cascades","rumors"),fill=c("red","yellow","blue"))
    dev.off()

    png("lambda_c_distrib_3D.png",width=2*480,height=400)
    par(mfrow=c(1,3))
    prior=allparameters.dataframe$lambda_c
    rg=range(prior)
    prior=density(prior,from=rg[1],to=rg[2])
    cols=rainbow(length(idscores))
    sapply(metrics,function(m){
           alldens=lapply(idscores,function(s)  density(allparameters.dataframe$lambda_c[rank(allscores[[s]][[m]],ties="first")<100],from=rg[1],to=rg[2]))
           yrange=range(sapply(alldens,function(d)range(d$y)),prior$y)
           plot(prior,xlab=expression(lambda_c),ylim=yrange,main=paste("Posterior lambda score=",m),log="x")
           lapply(1:length(idscores),function(i)lines(alldens[[i]],col=cols[i],lwd=3))
           legend("bottom",legend=idscores,col=cols,lwd=3)
    })
    dev.off()

    plot(density(posteriors[["kc"]][["size"]][["stime"]]))
    lines(density(allparameters.dataframe$lambda_c[rank(allscores[["kc"]][["size"]],ties="first")<500]))

for(param in colnames(allparameters.dataframe)){
    png(paste0(param,"_distrib_3D.png"),width=2*480,height=400)
    par(mfrow=c(1,3))
    prior=allparameters.dataframe[[param]]
    rg=range(prior)
    prior=density(prior,from=rg[1],to=rg[2])
    cols=rainbow(length(idscores))
    sapply(metrics,function(m){
           alldens=lapply(idscores,function(s)  density(posteriors[[s]][[m]][[param]],from=rg[1],to=rg[2]))
           yrange=range(sapply(alldens,function(d)range(d$y)),prior$y)
           plot(prior,xlab=param,ylim=yrange,main=paste("Posterior",param," score=",m))
           lapply(1:length(idscores),function(i)lines(alldens[[i]],col=cols[i],lwd=3))
           legend("bottom",legend=idscores,col=cols,lwd=3)
    })
    dev.off()
}
for(param in c("stime","dtime")){
    png(paste0(param,"_distrib_3D.png"),width=2*480,height=400)
    par(mfrow=c(1,3))
    prior=allparameters.dataframe[[param]]
    rg=range(prior)
    prior=density(prior,from=rg[1],to=rg[2],bw=1)
    cols=rainbow(length(idscores))
    sapply(metrics,function(m){
           alldens=lapply(idscores,function(s)  density(posteriors[[s]][[m]][[param]],from=rg[1],to=rg[2],bw=1))
           yrange=range(sapply(alldens,function(d)range(d$y)),prior$y)
           plot(prior,xlab=param,ylim=yrange,main=paste("Posterior",param," score=",m))
           lapply(1:length(idscores),function(i)lines(alldens[[i]],col=cols[i],lwd=3))
           legend("bottom",legend=idscores,col=cols,lwd=3)
    })
    dev.off()
}

    lapply(idscores,function(s)
           {
               lapply(mes,function(m)
                      {
                          png(paste0("param_corel_m=",m,"_s=",s,"_3D.png"),height=4*480,width=4*480)
                          plot(allparameters.dataframe[rank(allscores[[s]][[m]])<100,])
                          dev.off()
                      })
           })
    ###POSTCHECK FOR 3D
    prior_beta=c(-100,-10,0,10,100)
    prior_utility=c(-1,0,1)
    util=grep("utility.*",colnames(allparameters.dataframe))
    betad=grep("betaDistrib.*",colnames(allparameters.dataframe))
    s="kc"
    png("rumors_best_postcheck_3D.png")
    par(mfrow=c(1,3))
    sapply(metrics,function(m){
           plotCCFD(allca[[m]],main=m)
           best_param=unlist(allparameters.dataframe[which(rank(allscores[[s]][[m]]) == 1),])
           na=replicate(100,pointsCCFD(
                                       cascades3D(
                                                  N=best_param["N"],
                                                  stime=best_param["stime"],
                                                  dtime=best_param["dtime"],
                                                  IC=best_param["IC"],
                                                  Nmax=best_param["Nmax"],
                                                  R=best_param["R"],
                                                  time=best_param["repetition"],
                                                  captl=best_param["captl"],
                                                  utility = generateDistribeFromFrequencies(prior_utility,best_param["R"],best_param[util]),
                                                  betadistrib = generateDistribeFromFrequencies(prior_beta,best_param["N"],best_param[betad])
                                                  )[[m]],type="l",col=alpha("green",.1)))
           legend("bottomleft",legend=c("data","simulation using best fit"),col=c("black",alpha("red",.1)),lty=c(1,1),lwd=c(3,1))
           })
    dev.off()

    png("rumors_best_postcheck_3D_sub0_wrapper.png")
    par(mfrow=c(1,3))
    sapply(metrics,function(m){
           plotCCFD(allca[[m]],main=m)
           best_param=unlist(allparameters.dataframe[which(allscores[[s]][[m]] == min(allscores[[s]][[m]][allscores[[s]][[m]] > 0],na.rm=T)),])
           na=replicate(100,pointsCCFD(cascades3D.list(best_param )[[m]],type="l",col=alpha("green",.1)))
           legend("bottomleft",legend=c("data","simulation using best fit"),col=c("black",alpha("red",.1)),lty=c(1,1),lwd=c(3,1))
           })
    dev.off()
                                       test=cascades3D(
                                                  N=best_param["N"],
                                                  stime=best_param["stime"],
                                                  dtime=best_param["dtime"],
                                                  IC=best_param["IC"],
                                                  Nmax=best_param["Nmax"],
                                                  R=best_param["R"],
                                                  time=best_param["repetition"],
                                                  captl=best_param["captl"],
                                                  utility = generateDistribeFromFrequencies(prior_utility,best_param["R"],best_param[util]),
                                                  betadistrib = generateDistribeFromFrequencies(prior_beta,best_param["N"],best_param[betad])
                                                  )

    png("rumors_postcheck_3D_sub0.png")
    par(mfrow=c(1,3))
    sapply(metrics,function(m){
           plotCCFD(allca[[m]],main=m)
           apply(allparameters.dataframe,1,function(best_param){
                 best_param=unlist(best_param)
                 na=rep(5,pointsCCFD(
                                       cascades3D(
                                                  N=best_param["N"],
                                                  stime=best_param["stime"],
                                                  dtime=best_param["dtime"],
                                                  IC=best_param["IC"],
                                                  Nmax=best_param["Nmax"],
                                                  R=best_param["R"],
                                                  time=best_param["repetition"],
                                                  captl=best_param["captl"],
                                                  utility = generateDistribeFromFrequencies(prior_utility,best_param["R"],best_param[util]),
                                                  betadistrib = generateDistribeFromFrequencies(prior_beta,best_param["N"],best_param[betad])
                                                  )[[m]],type="l",col=alpha("green",.1)))
                 legend("bottomleft",legend=c("data","simulation using best fit"),col=c("black",alpha("red",.1)),lty=c(1,1),lwd=c(3,1))
                      })
           })
    dev.off()


}


abcNew <- function(){
    load("../results/testFull2/fullsim/res_sim_48.bin")
    load("../results/testFull1/scores.bin")
    load("../results/testFull1/parameters.bin")
    load("results/testFull1/scores.bin")
    load("../abcdir/observations.bin")
    plotCCFD(rud$rd$breadth)
    allfolders= list.dirs("newPrior/",recursive=F)

    mes=c("depth","breadth","size")
    kc=lapply(allfolders,function(j){load(paste0(j,"/scores.bin"));as.data.frame(t(sapply(scores,function(s)unlist(s[["kc"]][mes]))))})
    kc=do.call("rbind",kc)

    qc=lapply(allfolders,function(j)tryCatch({load(paste0(j,"/scores.bin"));as.data.frame(t(sapply(scores,function(s)unlist(s[["qc"]][mes]))))},error=function(e)return(NA)))
    qc=do.call("rbind",qc)
    kc=lapply(allfolders,function(j)tryCatch({load(paste0(j,"/scores.bin"));as.data.frame(t(sapply(scores,function(s)unlist(s[["kc"]][mes]))))},error=function(e)return(NA)))

    param=lapply(allfolders,function(j)as.data.frame({load(paste0(j,"/parameters.bin"));parameters}))
    parameters.df=do.call("rbind",param) 

    posteriors.kc=lapply(metrics,function(m,r=500) parameters.df[rank(kc[,m]) < r,])
    names(posteriors.kc)=metrics
    posteriors.qc=lapply(metrics,function(m,r=500) parameters.df[rank(qc[,m]) < r,])
    names(posteriors.qc)=metrics

par(mfrow=c(4,5),mar=c(5,5,1,1))
for(param in colnames(parameters.df)){
    prior=density(parameters.df[,param]) #prior distribution
    postC=density(posteriors.kc[["breadth"]][,param])  #posterior distribution for distance to cascade size
    rangex=range(parameters.df[,param])
    rangey=range(prior$y,postC$y)
    plot(prior,ylim=rangey,xlim=rangex,lwd=2,main=paste("Prior vs Posterior for param",param),xlab=param)
    lines(postC,col="yellow",lwd=2)
}

param="stime"
prior=density(parameters.df[,param],bw=1) #prior distribution
postC=density(posteriors.qc[["size"]][,param],bw=1)  #posterior distribution for distance to cascade size
rangex=range(parameters.df[,param])
rangey=range(prior$y,postC$y)
plot(prior,ylim=rangey,xlim=rangex,lwd=2,main=paste("Prior vs Posterior for param",param),xlab=param)
lines(postC,col="yellow",lwd=2)

par(mfrow=c(1,2))
median_best_param=apply(posteriors.kc$size,2,median)
plotCCFD(allca$size,main="median",xlab="cascades size")
na=replicate(10,pointsCCFD(
                            cascades3D(median_best_param)$size,type="l",col=alpha("red",.1)))
legend("bottomleft",legend=c("data","simulation ran median of posterior"),col=c("black",alpha("red",.1)),lty=c(1,1),lwd=c(3,1))

cols=topo.colors(length(seq(0.001,0.01,0.1)))

best_param=unlist(parameters.df[which(rank(qcrscors$size) == 1),])
plotCCFD(allca$size,main="best",xlab="cascades size")
na=replicate(10,pointsCCFD(cascades3D.list(best_param)$size,type="l",col=alpha("red",.1)))
legend("bottomleft",legend=c("data","simulation using best fit"),col=c("black",alpha("red",.1)),lty=c(1,1),lwd=c(3,1))


}

abcJunk <- function(){
    rud=c()
    print(rud)
    {
        newenv=new.env()
        allscores=lapply(1:2500,function(i){tryCatch({load(paste0("simulationbis/lastbsc_res_sim_",i,".bin"));return(list(score=rud$score,sum=rud$rd,time=rud$t))},error=function(err)return(NA))})
        print(rud)
        notlast=lapply(1:2500,function(i){tryCatch({load(paste0("simulationbis/notlastbsc_res_sim_",i,".bin"));return(list(score=rud$score,sum=rud$rd,time=rud$t))},error=function(err)return(NA))})
    }
    print(rud)
    max(unlist(alltimes),na.rm=T)
    score=as.numeric(unlist(sapply(allscores,function(alt)alt["score"])))
    times=lapply(allscores,function(alal){
                 if(!is.na(alal['time'])){
                     z=alal$time
                     units(z)="secs"
                     alal$time=z
                 }
})
    times=unlist(times)
    ulu=as.numeric(unlist(lapply(allscores,function(al)al["score"])))
    ulu=lapply(allscores,function(al)al$score)
    scores=order(unlist(ulu),na.last=F)
    plot(unlist(ulu)[scores])
    pointsCCFD(allcasca$size)
    pointsCCFD(allscores[[]]$sum$size,col="red")

    minimum=allscores[unlist(ulu)==min(ulu,na.rm=T) & !is.na(unlist(ulu))]
    maximum=allscores[unlist(ulu)==max(ulu,na.rm=T) & !is.na(unlist(ulu))]

    minimum=allscores[unlist(allnewsc)==min(allnewsc,na.rm=T) & !is.na(unlist(allnewsc))]
    maximum=allscores[unlist(allnewsc)==max(allnewsc,na.rm=T) & !is.na(unlist(allnewsc))]

    minimum=allscores[unlist(allKSscore)==min(allKSscore,na.rm=T) & !is.na(unlist(allKSscore))]
    maximum=allscores[unlist(allKSscore)==max(allKSscore,na.rm=T) & !is.na(unlist(allKSscore))]

    getMin <- function(listofexp,vecofscore) listofexp[unlist(vecofscore)==min(vecofscore,na.rm=T) & !is.na(unlist(vecofscore))]
    getMax <- function(listofexp,vecofscore) listofexp[unlist(vecofscore)==max(vecofscore,na.rm=T) & !is.na(unlist(vecofscore))]

    getMin <- function(listofexp) listofexp[unlist(vecofscore)==min(vecofscore,na.rm=T) & !is.na(unlist(vecofscore))]
    getMin <- function(listofexp) listofexp[lapply(scores,function(s)s$score)==min(sapply(scores,function(s)s$score))]
    getMin <- function(listofexp) listofexp[lapply(listofexp,function(s)s$score)==min(sapply(listofexp,function(s)s$score))]

    maximum=getMin(allscores,allKSscore)
    minimum=
    plotCCFD(allcasca$size,col=alpha(1,.2),type="o",cex=1)
    pointsCCFD(getMin(allscores,allKSscore)[[1]]$sum$size,col="green")
    pointsCCFD(,col="green")
    testuntruc=getMin(allscores,allQSscore)[[1]]$sum
    plotCCFD(tapply(testuntruc$size,testuntruc$rumor,sum))


    allKSscore=sapply(1:2500,function(i){tryCatch({load(paste0("simulationbis/lastbsc_res_sim_",i,".bin"));return(ks.test(rud$rd$size,allcasca$size)$statistic)},error=function(err)return(NA))})
    allQSscore=sapply(1:2500,function(i){tryCatch({load(paste0("simulationbis/lastbsc_res_sim_",i,".bin"));return(quantilediff(rud$rd$size,allcasca$size))},error=function(err)return(NA))})
    allQSscoreB=sapply(1:2500,function(i){tryCatch({load(paste0("simulationbis/notlastbsc_res_sim_",i,".bin"));return(quantilediff(rud$rd$size,allcasca$size))},error=function(err)return(NA))})
    allNscore=sapply(1:2500,function(i){tryCatch({load(paste0("simulationbis/lastbsc_res_sim_",i,".bin"));return(mean(diff(rud$rd$size,allcasca$size)))},error=function(err)return(NA))})
    allSumNscore=sapply(1:2500,function(i){tryCatch({load(paste0("simulationbis/lastbsc_res_sim_",i,".bin"));return(sum(diff(rud$rd$size,allcasca$size)))},error=function(err)return(NA))})
g

}

testalarache <- function(){
    rumors.scores$qc$size[1:2]
    repetead=apply(rumors.randposteriors[1:3,],1,function(param)replicate(10,quantilediff(getRumors(cascades3D.list(unlist(param))),allru$size)))

    getScore.list <- function(listp)
        quantilediff(getRumors(cascades3D.list(unlist(listp))),allru$size)

    randsample=sample.int(length(rumors.scores$qc$size),100)
    rumors.randposteriors=rumors.parameters[randsample,]
    rumors.randposteriors.scores=rumors.scores$qc$size[randsample]

    replicascade <- function(i,rep)
        cbind(rumors.scores$qc$size[i],replicate(rep,getScore.list(rumors.parameters[i,])))



repettion=do.call("rbind",lapply(sample(length(rumors.scores$qc$size),20),function(j)replicascade(j,100)))

best_param=unlist(rumors.parameters[which(rank(rumors.scores$qc$size) == 1),])
plotCCFD(allru$size,main="best",xlab="rumors size")
na=replicate(10,pointsCCFD(getRumors(cascades3D.list(best_param)),col=alpha("red",.1)))
legend("bottomleft",legend=c("data","simulation using best fit"),col=c("black",alpha("red",.1)),lty=c(1,1),lwd=c(3,1))



plotCCFD(allru$size,main="best",xlab="rumors size")
rrs=range(rumors.posteriors.scores.kc)
spacescore=round(seq(rrs[1],rrs[2],0.0000001),digit=7)
cols=topo.colors(length(spacescore))
names(cols)=spacescore
apply(rumors.posteriors$kc[order(rumors.posteriors.scores.kc,decreasing=T),],1,function(p)
      {
      best_param=p
      print(p)
      na=replicate(5,
                   {
                       rums=getRumors(cascades3D.list(unlist(p)))
                       sc=round(kldiff(rums,allru$size),digits = 7)
                       print(sc)
                       pointsCCFD(rums,col=alpha(cols[as.character(sc)],.6))
                   }
      )
      })
plotCCFD(allru$size,main="best",xlab="rumors size")
rrs=range(rumors.posteriors.scores.qc)
spacescore=round(seq(rrs[1],rrs[2],0.01),digit=2)
cols=topo.colors(length(spacescore))
names(cols)=spacescore
}


fromFolderToTest <- function(){

    

  plot(sort(vrefo.posteriors$qc$true$mu),ecdf(vrefo.posteriors$qc$true$mu)(sort(vrefo.posteriors$qc$true$mu)),col="green",type="l",log="x",xlab="mu",ylab="proba",main="CDF for posterior of mu")
  lines(sort(vrefo.posteriors$qc$false$mu),ecdf(vrefo.posteriors$qc$false$mu)(sort(vrefo.posteriors$qc$false$mu)),col="red")
     legend("topleft",legend=c("false","true"),col=c("red","green"),lwd=c(3,3))


    #With all the folder we can now get all the scores
    foldrum="../testneutralRumors"
    rumors.scores = getAllscores(foldrum,idscores = c("qc"),metrics=c("false","true","size"),lim=0:290)
    rumors.parameters = getAllparameters(foldrum,lim=0:290)
    rumors.posteriors=getAllposteriors(rumors.scores,rumors.parameters,500)
    rumors.posteriors.scores.qc=lapply(rumors.scores$qc,function(i)i[rank(i,ties="first")<500])


    foldrum="../testAllvsTFsinC/"
    sinC.scores = getAllscores(foldrum,idscores = c("qc"),lim=checkIfNumberGood(foldrum),metrics=c("size","true","false"),log=T)
    sinC.scores = getAllscores(foldrum,idscores = c("qc"),lim=0:6,metrics=c("size","true","false"),log=T)
    sinC.parameters = getAllparameters(foldrum,lim=checkIfNumberGood(foldrum),log=T)
    sinC.parameters = getAllparameters(foldrum,lim=0:6,log=T)
    #sinC.scores = getAllscores(foldrum,idscores = c("qc"),lim=(0:1149),metrics=c("size","true","false"))
    #sinC.parameters = getAllparameters(foldrum,lim=(0:1149))
    sinC.posteriors=getAllposteriors(sinC.scores,sinC.parameters,100)
    sinC.posteriors.scores.qc=lapply(sinC.scores$qc,function(i)i[rank(i,ties="first")<1000])

    foldrum="../testAllvsTF/"
    vrefo.scores = getAllscores(foldrum,idscores = c("qc"),lim=checkIfNumberGood(foldrum),metrics=c("true","false"),log=T)
    vrefo.parameters = getAllparameters(foldrum,lim=checkIfNumberGood(foldrum),log=T)
    vrefo.posteriors=getAllposteriors(vrefo.scores,vrefo.parameters,1000)
    vrefo.posteriors.scores.qc=lapply(vrefo.scores$qc,function(i)i[rank(i,ties="first")<1000])


    foldrum="../testAllvsTFMine"
    mine.scores = getAllscores(foldrum,lim=0:2,idscores = c("qc"),metrics=c("size","true","false"))
    mine.parameters = getAllparameters(foldrum,lim=0:2)
    mine.posteriors=getAllposteriors(mine.scores,mine.parameters,1000)
    mine.posteriors.scores.qc=lapply(mine.scores$qc,function(i)i[rank(i,ties="first")<1000])

    foldrum="../testAllvsTFMine2"
    mine.scores = getAllscores(foldrum,idscores = c("qc"),metrics=c("size","true","false"))
    mine.parameters = getAllparameters(foldrum)
    mine.posteriors=getAllposteriors(mine.scores,mine.parameters,1000)
    mine.posteriors.scores.qc=lapply(mine.scores$qc,function(i)i[rank(i,ties="first")<1000])

    foldrum="../testAllvsTFAlberto2"
    alberto.scores = getAllscores(foldrum,idscores = c("qc"),metrics=c("size"))
    alberto.parameters = getAllparameters(foldrum)
    alberto.posteriors=getAllposteriors(alberto.scores,alberto.parameters,1000)
    alberto.posteriors.scores.qc=lapply(alberto.scores$qc,function(i)i[rank(i,ties="first")<1000])

    foldrum="../testAllNeutral/"
    neutral.scores = getAllscores(foldrum,idscores = c("qc"),metrics=c("size"))
    neutral.parameters = getAllparameters(foldrum)
    neutral.posteriors=getAllposteriors(neutral.scores,neutral.parameters,500)
    neutral.posteriors.scores.qc=lapply(neutral.scores$qc,function(i)i[rank(i,ties="first")<500])

    foldrum="../testneutralRumorsVF/"
    lim=checkIfNumberGood(foldrum)
    neutral.scores = getAllscores(foldrum,idscores = c("qc"),lim=lim,metrics=c("true","false"))
    neutral.parameters = getAllparameters(foldrum,lim=lim)
    neutral.posteriors=getAllposteriors(neutral.scores,neutral.parameters,1000)
    neutral.posteriors.scores.qc=lapply(neutral.scores$qc,function(i)i[rank(i,ties="first")<1000])

    foldrum="../nonneutralRumorsVF/"
    lim=checkIfNumberGood(foldrum)
    nonneutral.scores = getAllscores(foldrum,idscores = c("qc"),lim=lim,metrics=c("true","false"))
    nonneutral.parameters = getAllparameters(foldrum,lim=lim)
    nonneutral.posteriors=getAllposteriors(nonneutral.scores,nonneutral.parameters,1000)
    nonneutral.posteriors.scores.qc=lapply(nonneutral.scores$qc,function(i)i[rank(i,ties="first")<1000])


    foldrum="../testAllvsTFAlbertoCorrected/"
    lim=checkIfNumberGood(foldrum)
    albertoc.scores = getAllscores(foldrum,idscores = c("qc"),lim=4:22,metrics=c("size","true","false"))
    albertoc.parameters = getAllparameters(foldrum,log=T,lim=5:22)
    albertoc.posteriors=getAllposteriors(albertoc.scores,albertoc.parameters,500)
    albertoc.posteriors.scores.qc=lapply(albertoc.scores$qc,function(i)i[rank(i,ties="first")<500])

    foldrum="../testAllvsTFAlberto"
    #rem=checkIfNumberGood(foldrum)
    #alberto.scores = getAllscores(foldrum,idscores = c("qc"),lim=(0:899)[-rem],metrics=c("size","true","false"))
    alberto.scores = getAllscores(foldrum,idscores = c("qc"),lim=0:400,metrics=c("size","true","false"))
    alberto.parameters = getAllparameters(foldrum,lim=0:400,log=T)
    #alberto.scores = getAllscores(foldrum,idscores = c("qc"),lim=(0:1149),metrics=c("size","true","false"))
    #alberto.parameters = getAllparameters(foldrum,lim=(0:1149))
    alberto.posteriors=getAllposteriors(alberto.scores,alberto.parameters,500)
    alberto.posteriors.scores.qc=lapply(alberto.scores$qc,function(i)i[rank(i,ties="first")<500])

    allexpes.scores=list(alberto=alberto.scores,mine=mine.scores,beta=vrefo.scores,random=sinC.scores)
    allexpes.parameters=list(alberto=alberto.parameters,mine=mine.parameters,beta=vrefo.parameters,random=sinC.parameters)

    nmin=min(sapply(allexpes.parameters,function(u)length(u[[1]]))) 

    lapply(allexpes.scores,function(u)lapply(u,function(v)lapply(v,function(w){
                                                                 print(paste("nmin: ",nmin))
                                                                 print(paste("length: ",length(w)))
                                                                 subsample(w,nmin)})))

    allexpes.scores.compare=lapply(allexpes.scores,lapply,lapply,subsample,nmin)
    allexpes.parameters.compare=lapply(allexpes.parameters,lapply,subsample,nmin)
    allexpes.posterios=lapply(names(allexpes.parameters.compare),function(n)getAllposteriors(allexpes.scores[[n]],allexpes.parameters[[n]],1000))
names(allexpes.posterios)=names(allexpes.parameters.compare)
    ## Plot posteriors 
    png("posteriors_param_random.png",width=1*480,height=1.5*480,pointsize = 16)

    cols=c(alpha("grey",.8),alpha("red",.8),alpha("green",.8))
    par(mfrow=c(3,2))
    plot2dens(A=alberto.posteriors$qc$false$mu,B=alberto.posteriors$qc$true$mu,prior=alberto.parameters$mu,xlab="mu",main="",from=0,to=0.001,cols=cols)
    legend("topright",legend=c("prior","false","true"),fill=cols,cex=1)
    plot2dens(A=alberto.posteriors$qc$false$Nmin,B=alberto.posteriors$qc$true$Nmin,prior=alberto.parameters$Nmin,xlab="Nmin",main="",cols=cols)
    plot2dens(A=alberto.posteriors$qc$false$t_step,B=alberto.posteriors$qc$true$t_step,prior=alberto.parameters$t_step,xlab="t_step",main="",cols=cols)
    plot2dens(A=alberto.posteriors$qc$false$tau,B=alberto.posteriors$qc$true$tau,prior=alberto.parameters$tau,xlab="tau",main="",cols=cols)
    plot2dens(A=alberto.posteriors$qc$false$C,B=alberto.posteriors$qc$true$C,prior=alberto.parameters$C,xlab="C",main="",cols=cols)
    plot2dens(A=alberto.posteriors$qc$false$TF,B=alberto.posteriors$qc$true$TF,prior=alberto.parameters$TF,xlab="tf",main="",cols=cols)
    dev.off()


    fromTestScoresVSParam(alberto.parameters,alberto.seores$qc$true,modelwrapper = randomCascades.list,repet = 50,numberparam = 4,data=trueru$size)
    plot(0,0,xlim=c(0,50),ylim=c(0,50),type="n",xlab="rerun",ylab="origin")
sapply(sample.int(nrow(alberto.posteriors$qc$true),500),function(i)replicate(10, points(quantilediff(randomCascades.list(unlist(alberto.posteriors$qc$true[i,]),alberto=T)$size,trueru$size),alberto.posteriors.scores.qc$true[i])))
sapply(sample.int(nrow(mine.posteriors$qc$true),100),function(i)replicate(10, points(quantilediff(randomCascades.list(unlist(mine.posteriors$qc$true[i,]),alberto=F)$size,trueru$size),mine.posteriors.scores.qc$true[i])))
sapply(sample.int(nrow(mine.parameters),100),function(i)replicate(10, points(quantilediff(randomCascades.list(unlist(mine.parameters[i,]),alberto=F)$size,trueru$size),mine.scores$qc$true[i])))

    #test if score / param are proportional
    par(mfrow=c(1,2))
    fromTestScoresVSParam(sinC.parameters,sinC.scores$qc$true,modelwrapper = randomCascades.list,repet = 50,numberparam = 4,data=trueru$size)
    fromTestScoresVSParam(vrefo.parameters,vrefo.scores$qc$true,modelwrapper = randomCascades.list,repet = 50,numberparam = 4,data=trueru$size)
    fromTestScoresVSParam(rumors.parameters,rumors.scores$qc$true,modelwrapper = cascades3D.list,rumor=T,repet = 50,numberparam = 4,data=trueru$size)
    fromTestScoresVSParam(mine.parameters,mine.scores$qc$true,modelwrapper = randomCascades.list,rumor=F,repet = 100,numberparam = 10,data=trueru$size)
    fromTestScoresVSParam(alberto.parameters,alberto.scores$qc$true,modelwrapper = randomCascades.list,rumor=F,repet = 100,numberparam = 10,data=trueru$size)

    
    ## Plot posteriors 
a   png("posteriors_param_random.png",width=1*480,height=1.5*480,pointsize = 16)
    cols=c(alpha("grey",.8),alpha("red",.8),alpha("green",.8))
    par(mfrow=c(2,2))
    #plot2dens(A=sinC.posteriors$qc$false$beta,B=sinC.posteriors$qc$true$beta,prior=sinC.parameters$beta,xlab="beta",main="",cols=cols)
    plot2dens(A=sinC.posteriors$qc$false$mu,B=sinC.posteriors$qc$true$mu,prior=sinC.parameters$mu,xlab="mu",main="",from=0,to=0.003,cols=cols)
    legend("topright",legend=c("prior","false","true"),fill=cols,cex=1)
    plot2dens(A=sinC.posteriors$qc$false$Nmin,B=sinC.posteriors$qc$true$Nmin,prior=sinC.parameters$Nmin,xlab="Nmin",main="",cols=cols)
    plot2dens(A=sinC.posteriors$qc$false$t_step,B=sinC.posteriors$qc$true$t_step,prior=sinC.parameters$t_step,xlab="t_step",main="",cols=cols)
    plot2dens(A=sinC.posteriors$qc$false$tau,B=sinC.posteriors$qc$true$tau,prior=sinC.parameters$tau,xlab="tau",main="",cols=cols)
    dev.off()

    ## Plot posteriors 
    png("posteriors_param_random.png",width=1*480,height=1.5*480,pointsize = 16)
    cols=c(alpha("grey",.8),alpha("red",.8),alpha("green",.8))
    par(mfrow=c(3,2))
    #plot2dens(A=mine.posteriors$qc$false$beta,B=mine.posteriors$qc$true$beta,prior=mine.parameters$beta,xlab="beta",main="",cols=cols)
    plot2dens(A=mine.posteriors$qc$false$mu,B=mine.posteriors$qc$true$mu,prior=mine.parameters$mu,xlab="mu",main="",from=0,to=0.002,cols=cols)
    legend("topright",legend=c("prior","false","true"),fill=cols,cex=1)
    plot2dens(A=mine.posteriors$qc$false$Nmin,B=mine.posteriors$qc$true$Nmin,prior=mine.parameters$Nmin,xlab="Nmin",main="",cols=cols)
    plot2dens(A=mine.posteriors$qc$false$t_step,B=mine.posteriors$qc$true$t_step,prior=mine.parameters$t_step,xlab="t_step",main="",cols=cols)
    plot2dens(A=mine.posteriors$qc$false$tau,B=mine.posteriors$qc$true$tau,prior=mine.parameters$tau,xlab="tau",main="",cols=cols)
    plot2dens(A=mine.posteriors$qc$false$TF,B=mine.posteriors$qc$true$TF,prior=mine.parameters$TF,xlab="TF",main="",cols=cols)
    plot2dens(A=mine.posteriors$qc$false$C,B=mine.posteriors$qc$true$C,prior=mine.parameters$C,xlab="C",from=0,to=0.08,main="",cols=cols)
    dev.off()

    png("posteriors_param_randomwithbeta.png",width=1.8*480,height=1.5*480,pointsize = 26)
    cols=c(alpha("grey",.8),alpha("red",.8),alpha("green",.8))
    par(mfrow=c(2,3))
    par(mar=c(4,4,1,1))
    plot2dens(A=vrefo.posteriors$qc$false$mu,B=vrefo.posteriors$qc$true$mu,prior=vrefo.parameters$mu,xlab="mu",main="",from=0,to=0.003,cols=cols)
    legend("topright",legend=c("prior","false","true"),fill=cols,cex=1)
    plot2dens(A=vrefo.posteriors$qc$false$Nmin,B=vrefo.posteriors$qc$true$Nmin,prior=vrefo.parameters$Nmin,xlab="Nmin",main="",cols=cols)
    plot2dens(A=vrefo.posteriors$qc$false$beta,B=vrefo.posteriors$qc$true$beta,prior=vrefo.parameters$beta,xlab="beta",main="",cols=cols)
    plot2dens(A=vrefo.posteriors$qc$false$t_step,B=vrefo.posteriors$qc$true$t_step,prior=vrefo.parameters$t_step,xlab="t_step",main="",cols=cols)
    plot2dens(A=vrefo.posteriors$qc$false$tau,B=vrefo.posteriors$qc$true$tau,prior=vrefo.parameters$tau,xlab="tau",main="",cols=cols)
    dev.off()

    png("posteriors_param_randomwithbetaVSnobeta_false.png",width=1*480,height=1.5*480,pointsize = 16)
    cols=c(alpha("grey",.8),alpha("red",.8),alpha("green",.8))
    par(mfrow=c(2,2),main="")
    plot2dens(A=sinC.posteriors$qc$false$mu,B=vrefo.posteriors$qc$false$mu,prior=sinC.parameters$mu,xlab="mu",main="false tweets",from=0,to=0.001,cols=cols)
    legend("topright",legend=c("prior","false no beta","false beta"),fill=cols,cex=1)
    plot2dens(A=sinC.posteriors$qc$false$Nmin,B=vrefo.posteriors$qc$false$Nmin,prior=sinC.parameters$Nmin,xlab="Nmin",main="false tweets",cols=cols)
    plot2dens(A=sinC.posteriors$qc$false$t_step,B=vrefo.posteriors$qc$false$t_step,prior=sinC.parameters$t_step,xlab="t_step",main="false tweets",cols=cols)
    plot2dens(A=sinC.posteriors$qc$false$tau,B=vrefo.posteriors$qc$false$tau,prior=sinC.parameters$tau,xlab="tau",main="false tweets",cols=cols)
    dev.off()
    
    png("posteriors_param_randomwithbetaVSnobeta_true.png",width=1*480,height=1.5*480,pointsize = 16)
    cols=c(alpha("grey",.8),alpha("red",.8),alpha("green",.8))
    par(mfrow=c(2,2),main="")
    plot2dens(A=sinC.posteriors$qc$true$mu,B=vrefo.posteriors$qc$true$mu,prior=sinC.parameters$mu,xlab="mu",main="true tweets",from=0,to=0.003,cols=cols)
    legend("topright",legend=c("prior","true no beta","true beta"),fill=cols,cex=1)
    plot2dens(A=sinC.posteriors$qc$true$Nmin,B=vrefo.posteriors$qc$true$Nmin,prior=sinC.parameters$Nmin,xlab="Nmin",main="true tweets",cols=cols)
    plot2dens(A=sinC.posteriors$qc$true$t_step,B=vrefo.posteriors$qc$true$t_step,prior=sinC.parameters$t_step,xlab="t_step",main="true tweets",cols=cols)
    plot2dens(A=sinC.posteriors$qc$true$tau,B=vrefo.posteriors$qc$true$tau,prior=sinC.parameters$tau,xlab="tau",main="true tweets",cols=cols)
    dev.off()

    png("posteriors_param_.png",width=1*480,height=1.5*480,pointsize = 16)
    cols=c(alpha("grey",.8),alpha("red",.8),alpha("green",.8))
    par(mfrow=c(2,2))
    plot2dens(A=sinC.posteriors$qc$true$mu,B=vrefo.posteriors$qc$true$mu,prior=sinC.parameters$mu,xlab="mu",main="true tweets",from=0,to=0.003,cols=cols)
    legend("topright",legend=c("prior","true no beta","true beta"),fill=cols,cex=1)
    plot2dens(A=sinC.posteriors$qc$true$Nmin,B=vrefo.posteriors$qc$true$Nmin,prior=sinC.parameters$Nmin,xlab="Nmin",main="true tweets",cols=cols)
    plot2dens(A=sinC.posteriors$qc$true$t_step,B=vrefo.posteriors$qc$true$t_step,prior=sinC.parameters$t_step,xlab="t_step",main="true tweets",cols=cols)
    plot2dens(A=sinC.posteriors$qc$true$tau,B=vrefo.posteriors$qc$true$tau,prior=sinC.parameters$tau,xlab="tau",main="true tweets",cols=cols)
    dev.off()
    ### plot reruns of posteriors 

    png("posteriors_rerurn_random.png",width=1.5*480)
    par(mfrow=c(1,2))
    allrz.sinC.false=plotPosteriorsCheck(sinC.posteriors$qc$false,sinC.posteriors.scores.qc$false,modelwrapper=randomCascades.list,data=falseru$size,scorefun=quantilediff,rep=50,samplepost=10000,clrs=alpha("red",.1),type="sampling",rumor=F,main="resampling from posterior: false tweet",xlab="size of cascade" )
    pointsCCFD(falseru$size,col=alpha("black",.3))
    allrz.sinC.true=plotPosteriorsCheck(sinC.posteriors$qc$true,sinC.posteriors.scores.qc$true,modelwrapper=randomCascades.list,data=trueru$size,scorefun=quantilediff,rep=50,samplepost=10000,clrs=alpha("green",.1),type="sampling",rumor=F,main="resampling from posterior: true tweet",xlab="size of cascade" )
    pointsCCFD(trueru$size,col=alpha("black",.3))
    dev.off()

    png("posteriors_rerurn_topfivemine.png",width=1.5*480)
    par(mfrow=c(1,2))
    allrz.mine.false=plotPosteriorsCheck(mine.posteriors$qc$false,mine.posteriors.scores.qc$false,modelwrapper=randomCascades.list,data=falseru$size,scorefun=quantilediff,rep=3,samplepost=10000,clrs=alpha("red",.2),type="sampling",rumor=F,main="resampling from posterior: false tweet",xlab="size of cascade" )
    pointsCCFD(falseru$size,col=alpha("black",.3))
    allrz.mine.true=plotPosteriorsCheck(mine.posteriors$qc$true,mine.posteriors.scores.qc$true,modelwrapper=randomCascades.list,data=trueru$size,scorefun=quantilediff,rep=3,samplepost=10000,clrs=alpha("green",.4),type="sampling",rumor=F,main="resampling from posterior: true tweet",xlab="size of cascade" )
    pointsCCFD(trueru$size,col=alpha("black",.3))
    dev.off()
    png("posteriors_rerurn_topfivealberto.png",width=1.5*480)
    par(mfrow=c(1,2))

    allrz.alberto.false=plotPosteriorsCheck(alberto.posteriors$qc$false,alberto.posteriors.scores.qc$false,modelwrapper=randomCascades.list,data=falseru$size,scorefun=quantilediff,rep=3,samplepost=1000,clrs=alpha("red",.2),type="sampling",rumor=F,main="resampling from posterior: false tweet",xlab="size of cascade",alberto=T)
    pointsCCFD(falseru$size,col=alpha("black",.3))
    allrz.alberto.true=plotPosteriorsCheck(alberto.posteriors$qc$true,alberto.posteriors.scores.qc$true,modelwrapper=randomCascades.list,data=trueru$size,scorefun=quantilediff,rep=3,samplepost=1000,clrs=alpha("green",.4),type="sampling",rumor=F,main="resampling from posterior: true tweet",xlab="size of cascade" ,alberto=T)
    pointsCCFD(trueru$size,col=alpha("black",.3))
    dev.off()

    allrz.alberto.false=reruns(alberto.posteriors$qc$false,alberto.posteriors.scores.qc$false,modelwrapper=randomCascades.list,data=falseru$size,scorefun=quantilediff,rep=1,samplepost=1000,type="sampling",alberto=T)
    allrz.alberto.true=reruns(alberto.posteriors$qc$true,alberto.posteriors.scores.qc$true,modelwrapper=randomCascades.list,data=trueru$size,scorefun=quantilediff,rep=1,samplepost=1000,type="sampling",alberto=T)
    allrz.mine.false=reruns(mine.posteriors$qc$false,mine.posteriors.scores.qc$false,modelwrapper=randomCascades.list,data=falseru$size,scorefun=quantilediff,rep=1,samplepost=1000,type="sampling",mine=T)
    allrz.mine.true=reruns(mine.posteriors$qc$true,mine.posteriors.scores.qc$true,modelwrapper=randomCascades.list,data=trueru$size,scorefun=quantilediff,rep=1,samplepost=1000,type="sampling",mine=T)


    png("posteriors_rerurn_randomwithbeta.png",width=1.5*480)
    par(mfrow=c(1,2))
    allrz.beta.false=plotPosteriorsCheck(vrefo.posteriors$qc$false,vrefo.posteriors.scores.qc$false,modelwrapper=randomCascades.list,data=falseru$size,scorefun=quantilediff,rep=50,samplepost=10000,clrs=alpha("red",.2),type="sampling",rumor=F,main="resampling from posterior: false tweet",xlab="size of cascade" )
    pointsCCFD(falseru$size,col=alpha("black",.3))
    allrz.beta.true=plotPosteriorsCheck(vrefo.posteriors$qc$true,vrefo.posteriors.scores.qc$true,modelwrapper=randomCascades.list,data=trueru$size,scorefun=quantilediff,rep=50,samplepost=10000,clrs=alpha("green",.4),type="sampling",rumor=F,main="resampling from posterior: true tweet",xlab="size of cascade" )
    pointsCCFD(trueru$size,col=alpha("black",.3))
    dev.off()

    png("posteriors_rerurn_randomwithbeta.png",width=1.5*480)
    par(mfrow=c(1,2))
    allrz.beta.false=plotPosteriorsCheck(vrefo.posteriors$qc$false,vrefo.posteriors.scores.qc$false,modelwrapper=randomCascades.list,data=falseru$size,scorefun=quantilediff,rep=50,samplepost=10000,clrs=alpha("red",.2),type="sampling",rumor=F,main="resampling from posterior: false tweet",xlab="size of cascade" )
    pointsCCFD(falseru$size,col=alpha("black",.3))
    allrz.beta.true=plotPosteriorsCheck(vrefo.posteriors$qc$true,vrefo.posteriors.scores.qc$true,modelwrapper=randomCascades.list,data=trueru$size,scorefun=quantilediff,rep=50,samplepost=10000,clrs=alpha("green",.4),type="sampling",rumor=F,main="resampling from posterior: true tweet",xlab="size of cascade" )
    pointsCCFD(trueru$size,col=alpha("black",.3))
    dev.off()

    png("posteriors_rerurn_random3D.png",width=1.5*480)
    par(mfrow=c(1,2))
    allrz=plotPosteriorsCheck(rumors.posteriors$qc$false,rumors.posteriors.scores.qc$false,modelwrapper=cascades3D.list,data=falseru$size,scorefun=quantilediff,rep=50,samplepost=100,clrs=alpha("red",.1),type="sampling",rumor=F,main="resampling from posterior: false tweet",xlab="size of cascade" )
    pointsCCFD(falseru$size,col=alpha("black",.3))
    allrz=plotPosteriorsCheck(rumors.posteriors$qc$true,rumors.posteriors.scores.qc$true,modelwrapper=cascades3D.list,data=trueru$size,scorefun=quantilediff,rep=50,samplepost=100,clrs=alpha("green",.1),type="sampling",rumor=T,main="resampling from posterior: true tweet",xlab="size of cascade" )
    pointsCCFD(trueru$size,col=alpha("black",.3))
    dev.off()



    ### rerun simulation from posteriors
    rerunscores.false=parApply(cl,vrefo.posteriors$qc$false[sample.int(nrow(vrefo.posteriors$qc$false),1000,replace = T),],1,function(p)quantilediff(randomCascades.list(unlist(p))$size,falseru$size))
    rerunscores.true=parApply(cl,vrefo.posteriors$qc$true[sample.int(nrow(vrefo.posteriors$qc$true),1000,replace = T),],1,function(p)quantilediff(randomCascades.list(unlist(p))$size,trueru$size))
    rerunscores.size=parApply(cl,vrefo.posteriors$qc$size[sample.int(nrow(vrefo.posteriors$qc$size),1000,replace = T),],1,function(p)quantilediff(randomCascades.list(unlist(p))$size,allru$size))

    #rerunscores.false.a=c(rerunscores.false,rerunscores.false.a)
    #rerunscores.true.a=c(rerunscores.true,rerunscores.true.a)
    #rerunscores.size.a=c(rerunscores.size,rerunscores.size.a)

    cl<-makeCluster(4,type="FORK")
    rerunscores.sinC.false=parApply(cl,sinC.posteriors$qc$false[sample.int(nrow(sinC.posteriors$qc$false),1000,replace = T),],1,function(p)quantilediff(randomCascades.list(unlist(p))$size,falseru$size))
    rerunscores.sinC.true=parApply(cl,sinC.posteriors$qc$true[sample.int(nrow(sinC.posteriors$qc$true),1000,replace = T),],1,function(p)quantilediff(randomCascades.list(unlist(p))$size,trueru$size))
    ### rerun simulation from posteriors
    rerunscores.beta.false=parApply(cl,vrefo.posteriors$qc$false[sample.int(nrow(vrefo.posteriors$qc$false),1000,replace = T),],1,function(p)quantilediff(randomCascades.list(unlist(p))$size,falseru$size))
    rerunscores.beta.true=parApply(cl,vrefo.posteriors$qc$true[sample.int(nrow(vrefo.posteriors$qc$true),1000,replace = T),],1,function(p)quantilediff(randomCascades.list(unlist(p))$size,trueru$size))
    stopCluster(cl)

    rerunscores.alberto.false=parApply(cl,alberto.posteriors$qc$false[sample.int(nrow(alberto.posteriors$qc$false),1000,replace = T),],1,function(p)quantilediff(randomCascades.list(unlist(p),alberto=T)$size,falseru$size))
    rerunscores.alberto.true=parApply(cl,alberto.posteriors$qc$true[sample.int(nrow(alberto.posteriors$qc$true),1000,replace = T),],1,function(p)quantilediff(randomCascades.list(unlist(p),alberto=T)$size,trueru$size))

    rerunscores.mine.false=parApply(cl,mine.posteriors$qc$false[sample.int(nrow(mine.posteriors$qc$false),1000,replace = T),],1,function(p)quantilediff(randomCascades.list(unlist(p),alberto=F)$size,falseru$size))
    rerunscores.mine.true=parApply(cl,mine.posteriors$qc$true[sample.int(nrow(mine.posteriors$qc$true),1000,replace = T),],1,function(p)quantilediff(randomCascades.list(unlist(p),alberto=F)$size,trueru$size))

    C.rerunscores.beta.false=sapply(allrz.bet.false,function(i)i[[1]]$score)
    C.rerunscores.beta.true=sapply(allrz.beta.true,function(i)i[[1]]$score)
    C.rerunscores.sinC.false=sapply(allrz.sinC.false,function(i)i[[1]]$score)
    C.rerunscores.sinC.true=sapply(allrz.sinC.true,function(i)i[[1]]$score)
    C.rerunscores.mine.false=sapply(allrz.mine.false,function(i)i[[1]]$score)
    C.rerunscores.mine.true=sapply(allrz.mine.true,function(i)i[[1]]$score)
    C.rerunscores.alberto.false=sapply(allrz.alberto.false,function(i)i[[1]]$score)
    C.rerunscores.alberto.true=sapply(allrz.alberto.true,function(i)i[[1]]$score)
    lapply(allrz,lapply,sapply(function(j)j$true$score))
    lapply(allrz,function(u)sapply(u,function(j)sapply(i,function(i)i[[1]]$score)))
    scores=lapply(allrz,function(exp)lapply(exp,function(tf)sapply(tf,function(u)u[[1]]$score)))

    rangex=range(sapply(scores,sapply,range))
    alldensities=lapply(scores,lapply,density)
    rangey=range(lapply(alldensities,lapply,"[[","y" ))
    explty=1:4
    names(explty)=names(alldensities)
    plot(1,1,ylim=rangey,xlim=rangex,type="n")
    for(n in names(alldensities))
        for(t in names(tf))
            lines(alldensities[[n]][[t]],col=tfcols[t],lty=explty[n])

     legend("toprigh",legend=c(tf,names(alldensities)),col=c("red","green",1,1,1,1),lwd=2,lty=c(1,1,1:4))

     plot(density(c(C.rerunscores.mine.true)),col="green",lwd=3,xlim=c(0,50),main="distance of fitted model to data",ylab="distance to data")

    png("compareBothModel.png")
     plot(density(c(C.rerunscores.beta.true,rerunscores.beta.true)),col="green",lwd=3,xlim=c(0,50),main="distance of fitted model to data",ylab="distance to data")
     lines(density(c(C.rerunscores.beta.false,rerunscores.beta.false)),col="red",lwd=3)
     lines(density(c(C.rerunscores.beta.true,rerunscores.beta.true)),col="green",lwd=3)
     lines(density(c(C.rerunscores.sinC.false,rerunscores.sinC.false)),col="red",lwd=3,lty=2)
     lines(density(c(C.rerunscores.sinC.true,rerunscores.sinC.true)),col="green",lwd=3,lty=2)
     lines(density(C.rerunscores.alberto.false),col="red",lwd=3,lty=3)
     lines(density(C.rerunscores.alberto.true),col="green",lwd=3,lty=3)
     lines(density(C.rerunscores.mine.false),col="red",lwd=3,lty=4)
     lines(density(C.rerunscores.mine.true),col="green",lwd=3,lty=4)
     legend("toprigh",legend=c("false","true","neutral","conformity"),col=c("red","green",1,1),lwd=2,lty=c(1,1,1,2))
     dev.off()

}
    

computeFrequenciesofAparition <- function(){
     ########COMPUTE THE PROBA OF APPARITION IN ONE DAY
     par(mfrow=c(1,2),mar=c(5,5,1,1))
     rateintro=sort(tapply(allca$days-min(allca$days),allca$rumor_id,min)) #for each rumor we assign the date of the first cascade that's spread it.
     alldayz=seq(min(rateintro),max(rateintro),3600*24) #generate all possible days, will be used to avoid skipping the day whith not new rumors in order to to take into account the day where number of rumor = 0
     freqD=as.numeric(table(factor(rateintro,levels=alldayz)))
     plot(as.POSIXct(alldayz,origin=min(allca$days),tz="GMT"),freqD,type="l",xlab="days",ylab="new rumors",cex=.8)     
     plot(density(freqD,bw=1,from=0),xlab="Number new rumors during one day",main="")


     alldate=as.POSIXct(trueca$date)
     trueca$hours=strptime(alldate,format="%Y-%m-%d %H", tz="GMT") #trunc by hours
     trueca$days=strptime(alldate,format="%Y-%m-%d", tz="GMT")     #trunc by days
     rateintro.true=sort(tapply(trueca$days-min(trueca$days),trueca$rumor_id,min)) #for each rumor we assign the date of the first cascade that's spread it.
     alldayz.true=seq(min(rateintro.true),max(rateintro.true),3600*24) #generate all possible days, will be used to avoid skipping the day whith not new rumors in order to to take into account the day where number of rumor = 0
     freqD.true=as.numeric(table(factor(rateintro.true,levels=alldayz.true)))

     alldate=as.POSIXct(falseca$date)
     falseca$hours=strptime(alldate,format="%Y-%m-%d %H", tz="GMT") #trunc by hours
     falseca$days=strptime(alldate,format="%Y-%m-%d", tz="GMT")     #trunc by days
     rateintro.false=sort(tapply(falseca$days-min(falseca$days),falseca$rumor_id,min)) #for each rumor we assign the date of the first cascade that's spread it.
     alldayz.false=seq(min(rateintro.false),max(rateintro.false),3600*24) #generate all possible days, will be used to avoid skipping the day whith not new rumors in order to to take into account the day where number of rumor = 0
     freqD.false=as.numeric(table(factor(rateintro.false,levels=alldayz.false)))

     plot(density(freqD.true,bw=1,from=0),xlab="Number new rumors during one day",main="")
     lines(density(freqD.false,bw=1,from=0),col="red",lwd=3)
     lines(density(freqD.true,bw=1,from=0),col="green",lwd=3)
     legend("bottomleft",legend=c("true","false"),col=c("red","green"),lwd=c(3,3))
}

###I guess all is in finalgraphs.R
binCountPlot <- function(){

    ###Plot binned distriubiton of real data (and not CCFD)
    plot(bins[1:(length(bins)-1)],as.numeric(prop.table(table(cut(allca$size,breaks=bins)))),log="xy",col=1,pch=20,type="l")
    points(bins[1:(length(bins)-1)],as.numeric(prop.table(table(cut(falseca$size,breaks=bins)))),col="red",pch=20)
    points(bins[1:(length(bins)-1)],as.numeric(prop.table(table(cut(trueca$size,breaks=bins)))),col="green",pch=20)

}




testconformity <- function(){
    plotCCFD(falseru$size)
    test=randomCascades(Nmin = 5000,t_steps = 100,mu = 0.00005, tau =25 , conformity=F,topfive = T,TF=5,C=.001,beta=0,alberto=F)

    par(mfrow=c(2,2))
    sapply(seq(2,10,length.out = 4),function(tf){
    plotCCFD(falseru$size)
    replicate(3,{
    test=randomCascades(Nmin = 5000,t_steps = 100,mu = 0.00005, tau =25 , conformity=F,topfive = T,TF=tf,C=.1,beta=0,alberto=F)
    pointsCCFD(test$size,col="grey")
      })
    replicate(3,{
    test=randomCascades(Nmin = 5000,t_steps = 100,mu = 0.00005, tau =25 , conformity=F,topfive = T,TF=tf,C=.1,beta=0,alberto=T)
    pointsCCFD(test$size,col="orange")
      })
    })

    test=randomCascades(Nmin = 5000,t_steps = 100,mu = 0.00005, tau =25 , conformity=F,topfive = T,TF=tf,C=.01,beta=0,alberto=T)

    microbenchmark(randomCascades(Nmin = 5000,t_steps = 100,mu = 0.00005, tau =25 , conformity=F,topfive = T,TF=5,C=.001,beta=0,alberto=T),randomCascades(Nmin = 5000,t_steps = 100,mu = 0.00005, tau =25 , conformity=F,topfive = T,TF=5,C=.001,beta=0,alberto=F),times = 10)

   pop= table(trueru)
   TF=5
   topfive=pop[rank(pop,ties.method = "first")>(length(pop)-TF)]
   topfive=as.numeric(names(topfive))
   C=.1
   cpop=sample(topfive,length(pop)*C,replace=T)
    
}


moreDraftTest <- function(){
    
    spl=sample(nrow(sinC.parameters),1000)

    simsN=lapply(spl,function(i)randomCascades.list(unlist(mine.parameters[i,])))
    simsN=lapply(spl,function(i)randomCascades.list(unlist(sinC.parameters[i,])))
    simsN=apply(sinC.posteriors$qc$true,1,function(i)randomCascades.list(i)))

    setdiff=as.data.frame(t(sapply(simsN,function(sn)c(ks.test(sn$size,trueru$size)$statistic,euclidediff(sn$size,trueru$size,nbin = 100),euclidediff(sn$size,trueru$size,nbin = 50),euclidediff(sn$size,trueru$size),quantilediff(sn$size,trueru$size)))))
    colnames(setdiff)=c("KS","BDiff100","BDiff50","BDiff20","QDiff")
    #pdf("compareCorrelationDistanceKSQDIFFBDIFF.pdf")
    plot(setdiff$KS ~ setdiff$QDiff,xlab="quantile diff", ylab="KS")
    #dev.off()

    png("compareBothModelAndDataTOPFIVE.png",width=2.3*480,height=2.2*480,pointsize = 20)
    par(mfrow=c(2,2))

    bins=10^(seq(0,6,.1))
    par(mar=c(3,3,1,1))
    allbined.false=lapply(allrz.mine.false,function(l)as.numeric(prop.table(table(cut(l[[1]]$distrib,breaks=bins,include.lowest = T)))))
    dfallbin.false=do.call("rbind",allbined.false)
    dfallbin.false[dfallbin.false==0]=0.0001
    dfallbin.false.log=log10(dfallbin.false)
    hdrcde::hdr.boxplot(lapply(1:ncol(dfallbin.false),function(i)dfallbin.false.log[,i]),space = .1,col=alpha("red",c(.1,.3)),ylim=c(-4,0),prob = c(50,95),yaxt="n",outline = T,pch=".")
    points(log10(as.numeric(prop.table(table(cut(falseru$size,breaks=bins,include.lowest = T))))),col="red",pch=20)
    points(log10(as.numeric(prop.table(table(cut(falseru$size,breaks=bins,include.lowest = T))))))
    axis(1,at = 1:length(bins),labels =  bins)
    axis(2,at = -4:0,labels =  c(0,10^(-3:0)))
    legend("bottomleft",lty=c(NA,NA,1,NA),pch=c(NA,NA,NA,21),fill=c(alpha("red",c(.1,.3)),NA,NA),border=c(1,1,NA,NA),legend=c("95% HDR","50% HDR","mode","data"),pt.bg=c(NA,NA,NA,"red"))

    par(mar=c(3,3,1,1))
    allbined.true=lapply(allrz.mine.true,function(l)as.numeric(prop.table(table(cut(l[[1]]$distrib,breaks=bins,include.lowest = T)))))
    dfallbin.true=do.call("rbind",allbined.true)
    dfallbin.true[dfallbin.true==0]=0.0001
    dfallbin.true.log=log10(dfallbin.true)
    hdrcde::hdr.boxplot(lapply(1:ncol(dfallbin.true),function(i)dfallbin.true.log[,i]),space = .1,col=alpha("green",c(.1,.3)),ylim=c(-4,0),prob = c(50,95),yaxt="n",outline = T,pch=".")
    points(log10(as.numeric(prop.table(table(cut(trueru$size,breaks=bins,include.lowest = T))))),col="green",pch=20)
    points(log10(as.numeric(prop.table(table(cut(trueru$size,breaks=bins,include.lowest = T))))))
    axis(1,at = 1:length(bins),labels =  bins)
    axis(2,at = -4:0,labels =  c(0,10^(-3:0)))
    legend("bottomleft",lty=c(NA,NA,1,NA),pch=c(NA,NA,NA,21),fill=c(alpha("green",c(.1,.3)),NA,NA),border=c(1,1,NA,NA),legend=c("95% HDR","50% HDR","mode","data"),pt.bg=c(NA,NA,NA,"green"))


    par(mar=c(3,3,1,1))
    allbined.false=lapply(allrz.alberto.false,function(l)as.numeric(prop.table(table(cut(l[[1]]$distrib,breaks=bins,include.lowest = T)))))
    dfallbin.false=do.call("rbind",allbined.false)
    dfallbin.false[dfallbin.false==0]=0.0001
    dfallbin.false.log=log10(dfallbin.false)
    hdrcde::hdr.boxplot(lapply(1:ncol(dfallbin.false),function(i)dfallbin.false.log[,i]),space = .1,col=alpha("red",c(.1,.3)),ylim=c(-4,0),prob = c(50,95),yaxt="n",outline = T,pch=".")
    points(log10(as.numeric(prop.table(table(cut(falseru$size,breaks=bins,include.lowest = T))))),col="red",pch=20)
    points(log10(as.numeric(prop.table(table(cut(falseru$size,breaks=bins,include.lowest = T))))))
    axis(1,at = 1:length(bins),labels =  bins)
    axis(2,at = -4:0,labels =  c(0,10^(-3:0)))
    legend("bottomleft",lty=c(NA,NA,1,NA),pch=c(NA,NA,NA,21),fill=c(alpha("red",c(.1,.3)),NA,NA),border=c(1,1,NA,NA),legend=c("95% HDR","50% HDR","mode","data"),pt.bg=c(NA,NA,NA,"red"))

    par(mar=c(3,3,1,1))
    allbined.true=lapply(allrz.alberto.true,function(l)as.numeric(prop.table(table(cut(l[[1]]$distrib,breaks=bins,include.lowest = T)))))
    dfallbin.true=do.call("rbind",allbined.true)
    dfallbin.true[dfallbin.true==0]=0.0001
    dfallbin.true.log=log10(dfallbin.true)
    hdrcde::hdr.boxplot(lapply(1:ncol(dfallbin.true),function(i)dfallbin.true.log[,i]),space = .1,col=alpha("green",c(.1,.3)),ylim=c(-4,0),prob = c(50,95),yaxt="n",outline = T,pch=".")
    points(log10(as.numeric(prop.table(table(cut(trueru$size,breaks=bins,include.lowest = T))))),col="green",pch=20)
    points(log10(as.numeric(prop.table(table(cut(trueru$size,breaks=bins,include.lowest = T))))))
    axis(1,at = 1:length(bins),labels =  bins)
    axis(2,at = -4:0,labels =  c(0,10^(-3:0)))
    legend("bottomleft",lty=c(NA,NA,1,NA),pch=c(NA,NA,NA,21),fill=c(alpha("green",c(.1,.3)),NA,NA),border=c(1,1,NA,NA),legend=c("95% HDR","50% HDR","mode","data"),pt.bg=c(NA,NA,NA,"green"))
    dev.off()


    abackupclearner <- function(){
        vrefo.posteriors
    }
    subsample <- function(x,y)x[sample.int(length(x),y)]

    plot(density(mine.posteriors$qc$false$mu* mine.posteriors$qc$false$t_step*mine.posteriors$qc$false$Nmin))
    plot(density(mine.posteriors$qc$true$mu* mine.posteriors$qc$true$t_step*mine.posteriors$qc$true$Nmin))



    listexp=c(
              "../testAllvsTFsinC/",
              "../testAllvsTF",
              "../testAllvsTFMine",
              "../testAllvsTFAlberto"
              )
    names(listexp)=c("random","conformist","topfiveS","topfiveA")

    il=lapply(listexp,getAllscores,idscores = c("qc"),metrics=c("true","false"))

    png("compareall.png",width=1000,height=1000)
    cols=1:4
    names(cols)=names(listexp)
    names(il)=names(listexp)
    par(mfrow=c(2,2))
    for(d in tf){
    plot(density(il[[1]]$qc[[d]]),xlim=c(0,5),ylim=c(0,.0005),type="n",main=d)
    lapply(names(il),function(u)lines(density(il[[u]]$qc[[d]][!is.na(il[[u]]$qc[[d]])]),col=cols[u],lwd=3))
    legend("topleft",col=cols,legend=names(cols),lwd=3) 
    }
    for(d in tf){
    plot(density(il[[1]]$qc[[d]]),xlim=c(0,50),ylim=c(0,.16),type="n",main=d)
    lapply(names(il),function(u)lines(density(il[[u]]$qc[[d]][!is.na(il[[u]]$qc[[d]])]),col=cols[u],lwd=3))
    legend("topleft",col=cols,legend=names(cols),lwd=3) 
    }
    dev.off()

    il.paramaters=lapply(listexp,getAllparameters,log=T)
    post=lapply(names(listexp),function(n)getAllposteriors(uuu[[n]],as.data.frame(uuu.param[[n]]),1000))
    names(post)=names(listexp[4])
    il.paramaters$topfiveS=lapply(il.paramaters$topfiveS,function(u)u[-which(is.na(il.paramaters$topfiveS$Nmin))])
    il$topfiveS$qc=lapply(il$topfiveS$qc,function(u)u[-which(is.na(u))])


   load(file="BACKUP")#il.paramaters,il)
   save(file="shorter_all",uuu.param,uuu)
   load(file="shorter_all")
    rm(il.paramaters,il)


    save(file="allrz_backB",allrz)


    size=list()
    size$true=trueru$size
    size$false=falseru$size
    bins=10^(seq(0,5.5,.1))
    tfcols=c("green","red")
    names(tfcols)=c("true","false")
    tf=c("true","false")
    names(tf)=c("true","false")

    rul=lapply(exp,function(u)u[unlist(lapply(u,lapply,function(i)if(is.na(i$score)) F else T))]) #remove NA 
   
    par(mfrow=c(1,2))
    llrz=list()
    bins=10^(seq(0,5.5,.1))
    lapply(allrz,function(exp){
           exp=allrz$topfiveA
        allbined=lapply(exp,lapply,function(l)as.numeric(prop.table(table(cut(l[[1]]$distrib,breaks=bins,include.lowest = T)))))
        na=lapply(tf,function(d)
                  {
                      allb=allbined[[d]]
                      par(mar=c(4,4,1,1))
                      dfallbin=do.call("rbind",allb)
                      dfallbin[dfallbin==0]=NA
                      dfallbin.log=log10(dfallbin)
                      hdrcde::hdr.boxplot(lapply((1:ncol(dfallbin))[-c(2,3,6)],function(i){
                                                 alog=dfallbin.log[,i]
                                                 alog[is.na(alog)]=-5
                                                 alog
                                                                 }
                      ),space = .1,col=alpha(tfcols[d],c(.1,.3)),ylim=c(-5.2,0),prob = c(50,95),yaxt="n",outline = T,lambda=1,h=.1,pch=".",xlab="Number of RTs",ylab="Frequencies")
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T)))))[-c(2,3,6)],col=tfcols[d],pch=20)
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T)))))[-c(2,3,6)])
                      axis(1,at = (1:length(bins))[0:6*10+1],labels =  bins[0:6*10+1])
                      axis(2,at = -5:0,labels =  c(0,10^(-4:0)))
                      legend(1,-3.5,lty=c(NA,NA,1,NA),pch=c(NA,NA,NA,21),fill=c(alpha(tfcols[d],c(.1,.3)),NA,NA),border=c(1,1,NA,NA),legend=c("95% HDR","50% HDR","mode","data"),pt.bg=c(NA,NA,NA,tfcols[d]))
                  }
        )
              }
        )

                      dfallbin=do.call("rbind",allbined$true)


    
plot2dens(,prior=log10(sample(1000:10000,100000,replace=T)/sample(100,100000,replace=T) ),from=-4,to=4)

plotNdens <- function(listdensA,lisdensB,listcols,listfill,from=NULL,to=NULL,prior=NULL,cols=c(alpha("red",.8),alpha("blue",.8),alpha("yellow",.8)),...){

    denseP=NULL
    denseA=NULL
    denseB=NULL
    if(!is.null(prior))prior=prior[!is.na(prior)]
    #if(!is.null(A))denseA=density(A,from=from,to=to)
    #if(!is.null(B))denseB=density(B,from=from,to=to)
    if(length(prior)==2)denseP=density(runif(5000,prior[1],prior[2]),from=from,to=to)
    else if(!is.null(prior))denseP=density(prior,from=from,to=to)

    for(t in tf){
        listdensA=lapply(post,function(u,v="mu",l=t)u$qc[[l]][[v]])
        to=.01
        from=min(sapply(listdensA,min),from)
        to=max(sapply(listdensA,max),to)

        alldensA=lapply(listdensA,function(a)density(a,from=from,to=0.2))
        listfill=1:4 
        names(listfill)=names(listdensA)
        names(cols)=c("P","A","B")
        rangex=range(sapply(alldensA,"[[","x"))
        rangey=range(sapply(alldensA,"[[","y"))
        plot(1,1,ylim=rangey,xlim=rangex,type="n")
        for(n in names(alldensA))
            polygon(c(from,alldensA[[n]]$x,to),c(0,alldensA[[n]]$y,0),col=alpha(tfcols[t],.5),lwd=2,density=listfill[n]*5,angle=45*listfill[n]/1.2)
        legend("topright",fill="red",density=listfill*5,angle=45*listfill/1.2,legend=names(alldensA),box.cex=c(2,1))
    }
        
    if(!is.null(A))
        polygon(c(from,denseA$x,to),c(0,denseA$y,0),col=cols["A"],lwd=2)
    if(!is.null(B))
        polygon(c(from,denseB$x,to),c(0,denseB$y,0),col=cols["B"],lwd=2)


}


checkResultsConsistantWithScore <- function(){
    allretest=list()
    for(exp in names(listexp)){
        h=sample.int(length(il.paramaters[[exp]]$Nmin),1)
        p=unlist(lapply(il.paramaters[[exp]],"[[",h))
        allretest[[exp]]=lapply(tf,function(l){
            s=il[[exp]]$qc[[l]][h]
            print(exp)
            alberto=F
            if(exp=="topfiveA")alberto=T
            s.tmp=replicate(50,quantilediff(randomCascades.list(p,alberto=alberto)$size,size[[l]]))
            return(list(os=s,ts=s.tmp))
        })
    }
    boxplot(allrete

            
            par(mfrow=c(2,2))
            for(exp in names(listexp))
                plot(t(sapply(sample(1000000,10),function(n){alberto=F;if(exp=="topfiveA")alberto=T;c(il[[exp]]$qc$true[n],quantilediff(randomCascades.list(unlist(lapply(il.paramaters[[exp]],"[[",n)),alberto=alberto)$size,size$true))})),main=exp)




}


addComputeRatios <- function(){
    cols=c(alpha("grey",.8),alpha("red",.8),alpha("green",.8))
    neutral.posteriors$qc$false$timestep =1.66*log(neutral.posteriors$qc$false$Nmax * neutral.posteriors$qc$false$repetition)-16
    neutral.posteriors$qc$true$timestep =1.66*log(neutral.posteriors$qc$true$Nmax * neutral.posteriors$qc$true$repetition)-16
    neutral.posteriors$qc$true$rat =log(neutral.posteriors$qc$true$N) / log(neutral.posteriors$qc$true$Nmax)
    neutral.posteriors$qc$false$rat =log(neutral.posteriors$qc$false$N) / log(neutral.posteriors$qc$false$Nmax)
    neutral.posteriors$qc$true$rrat =neutral.posteriors$qc$true$IC / neutral.posteriors$qc$true$R
    neutral.posteriors$qc$false$rrat =neutral.posteriors$qc$false$IC / neutral.posteriors$qc$false$R
    neutral.posteriors$qc$true$rcaptl =neutral.posteriors$qc$true$captl / neutral.posteriors$qc$true$IC
    neutral.posteriors$qc$false$rcaptl =neutral.posteriors$qc$false$captl / neutral.posteriors$qc$false$IC
    neutral.parameters$timestep =1.66*log(neutral.parameters$Nmax * neutral.parameters$repetition)-16
    neutral.parameters$rat =log(neutral.parameters$N) / log(neutral.parameters$Nmax)
    neutral.parameters$rrat =neutral.parameters$IC / neutral.parameters$R
    neutral.parameters$rcaptl =neutral.parameters$captl / neutral.parameters$IC
    par(mfrow=c(5,3))
    for(n in names(neutral.posteriors$qc$false))
        plot2dens(A=neutral.posteriors$qc$false[[n]],B=neutral.posteriors$qc$true[[n]],prior=neutral.parameters[[n]],xlab=n,main="",cols=cols)
    legend("topright",legend=c("prior","size","true"),fill=cols,cex=1)


    rerunbis[[exp]]=lapply(tf,function(l)reruns(cur[[l]],modelwrapper=cascades3D.list,data=size[[l]],scorefun=quantilediff,repet=0,samplepost=100,type="sampling",log=T))


}

##The following function are important to group the parameters with similiar meaning in order to have a visualisation that can be more easily understood 
expendConcatNonneutralProperties <- function(){
    util=grep("utility.*",colnames(nonneutral.parameters))
    betad=grep("betaDistrib.*",colnames(nonneutral.parameters))

    nonneutral.posteriors$qc$true$betaDistrib.new.neut = nonneutral.posteriors$qc$true$betaDistrib.new.2
    nonneutral.posteriors$qc$false$betaDistrib.new.neut = nonneutral.posteriors$qc$false$betaDistrib.new.2
    nonneutral.parameters$betaDistrib.new.neut = nonneutral.parameters$betaDistrib.new.2

    nonneutral.posteriors$qc$true$betaDistrib.new.pol = nonneutral.posteriors$qc$true$betaDistrib.new +nonneutral.posteriors$qc$true$betaDistrib.new.1 + nonneutral.posteriors$qc$true$betaDistrib.new.3 + nonneutral.posteriors$qc$true$betaDistrib.new.4 
    nonneutral.posteriors$qc$false$betaDistrib.new.pol = nonneutral.posteriors$qc$false$betaDistrib.new +nonneutral.posteriors$qc$false$betaDistrib.new.1 + nonneutral.posteriors$qc$false$betaDistrib.new.3 + nonneutral.posteriors$qc$false$betaDistrib.new.4 

    nonneutral.posteriors$qc$false$betaDistrib.new.neg = nonneutral.posteriors$qc$false$betaDistrib.new +nonneutral.posteriors$qc$false$betaDistrib.new.1 
    nonneutral.posteriors$qc$false$betaDistrib.new.pos = nonneutral.posteriors$qc$false$betaDistrib.new.3 + nonneutral.posteriors$qc$false$betaDistrib.new.4 

    nonneutral.posteriors$qc$true$betaDistrib.new.neg = nonneutral.posteriors$qc$true$betaDistrib.new +nonneutral.posteriors$qc$true$betaDistrib.new.1 
    nonneutral.posteriors$qc$true$betaDistrib.new.pos = nonneutral.posteriors$qc$true$betaDistrib.new.3 + nonneutral.posteriors$qc$true$betaDistrib.new.4 
    nonneutral.parameters$betaDistrib.new.neg = nonneutral.parameters$betaDistrib.new +nonneutral.parameters$betaDistrib.new.1 
    nonneutral.parameters$betaDistrib.new.pos = nonneutral.parameters$betaDistrib.new.3 +nonneutral.parameters$betaDistrib.new.4 

    nonneutral.parameters$betaDistrib.new.pol = nonneutral.parameters$betaDistrib.new +nonneutral.parameters$betaDistrib.new.1 + nonneutral.parameters$betaDistrib.new.3 + nonneutral.parameters$betaDistrib.new.4 

    nonneutral.posteriors$qc$true$utility.new.neut = nonneutral.posteriors$qc$true$utility.new.1
    nonneutral.posteriors$qc$false$utility.new.neut = nonneutral.posteriors$qc$false$utility.new.1
    nonneutral.parameters$utility.new.neut = nonneutral.parameters$utility.new.1

    nonneutral.posteriors$qc$true$utility.new.pol = nonneutral.posteriors$qc$true$utility.new + nonneutral.posteriors$qc$true$utility.new.2 
    nonneutral.posteriors$qc$false$utility.new.pol = nonneutral.posteriors$qc$false$utility.new + nonneutral.posteriors$qc$false$utility.new.2
    nonneutral.parameters$utility.new.pol = nonneutral.parameters$utility.new + nonneutral.parameters$utility.new.2

    nonneutral.posteriors$qc$false$utility.ratio = log(nonneutral.posteriors$qc$false$utility.new.pol /nonneutral.posteriors$qc$false$utility.new.neut)
    nonneutral.posteriors$qc$true$utility.ratio = log(nonneutral.posteriors$qc$true$utility.new.pol /nonneutral.posteriors$qc$true$utility.new.neut)
    nonneutral.parameters$utility.ratio = log(nonneutral.parameters$utility.new.pol /nonneutral.parameters$utility.new.neut)

    nonneutral.posteriors$qc$false$betaDistrib.ratio = log(nonneutral.posteriors$qc$false$betaDistrib.new.pol /nonneutral.posteriors$qc$false$betaDistrib.new.neut)
    nonneutral.posteriors$qc$true$betaDistrib.ratio = log(nonneutral.posteriors$qc$true$betaDistrib.new.pol /nonneutral.posteriors$qc$true$betaDistrib.new.neut)
    nonneutral.parameters$betaDistrib.ratio = log(nonneutral.parameters$betaDistrib.new.pol /nonneutral.parameters$betaDistrib.new.neut)
}


print2Dboxplot <- function(){
    #rerun[[exp]]=lapply(tf,function(l)reruns(cur[[l]],modelwrapper=cascades3D.list,data=size[[l]],scorefun=quantilediff,repet=0,samplepost=10,type="sampling",rumor=T,log=F))

    par(mfrow=c(3,4))
    dev.off()


    par(mfrow=c(2,2))
    par(mar=rep(1,4))

    prob=seq(50,100,10)
    hdr.boxplot.2d(nonneutral.parameters$betaDistrib.new.pol,nonneutral.parameters$utility.new.neut,prob=prob,xlim=c(0,1),ylim=c(0,1),shadecols=alpha("grey",.9),xlab="% polarized individual",ylab="% of neutral rumors",outside.points=F)
    for( i in rev(tf)){
        par(new=T)
        hdr.boxplot.2d(nonneutral.posteriors$qc[[i]]$betaDistrib.new.pol,nonneutral.posteriors$qc[[i]]$utility.new.neut,prob=prob,shadecols=tfcols[i],xlim=c(0,1),ylim=c(0,1),xlab="% polarized individual",ylab="% of neutral rumors",outside.points=F,show.points=F,pointcol=tfcols[i])
    }


    hdr.boxplot.2d(nonneutral.parameters$betaDistrib.new.neut,nonneutral.parameters$utility.new.neut,prob=prob,xlim=c(0,1),ylim=c(0,1),xlab="% neutral individual",ylab="% of neutral rumors",outside.points=F,pointcol=gfcols[i])
    for( i in (tf)){
        par(new=T)
        hdr.boxplot.2d(nonneutral.posteriors$qc[[i]]$betaDistrib.new.neut,nonneutral.posteriors$qc[[i]]$utility.new.neut,prob=prob,shadecols=tfcols[i],xlim=c(0,1),ylim=c(0,1),xlab="% neutral individual",ylab="% of neutral rumors",outside.points=F,pointcol=tfcols[i])
    }

    pdf("posteriorPoralizedRumorsVsAgents.pdf")
    hdr.boxplot.2d(nonneutral.parameters$utility.new.pol,nonneutral.parameters$betaDistrib.new.pol,prob=prob,shadecols=alpha("grey",.9),xlim=c(0,1),ylim=c(0,1)  ,ylab="% polarized individual",xlab="% of polarized rumors",outside.points=F,pointcol=tfcols[i])
    for( i in rev(tf)){
        par(new=T)
        hdrcde::hdr.boxplot.2d(nonneutral.posteriors$qc[[i]]$utility.new.pol,nonneutral.posteriors$qc[[i]]$betaDistrib.new.pol,prob=prob,shadecols=tfcols[i],xlim=c(0,1),ylim=c(0,1),ylab="% polarized individual",xlab="% of polarized rumors",outside.points=F,pointcol=tfcols[i])
    }
    dev.off()

    plot(density(log(nonneutral.parameters$utility.new.pol/nonneutral.parameters$betaDistrib.new.pol)),xlim=c(-2,2),ylim=c(0,2))
    for( i in rev(tf)){
        par(new=T)
        lines(density(log(nonneutral.posteriors$qc[[i]]$utility.new.pol/nonneutral.posteriors$qc[[i]]$betaDistrib.new.pol)),col=tfcols[i])
    }


    par(mfrow=c(1,2))
    hdr.boxplot.2d(nonneutral.parameters$betaDistrib.new.neg,nonneutral.parameters$utility.new.2,prob=prob,shadecols=alpha("grey",.9),xlim=c(0,1),ylim=c(0,1)  ,xlab="",ylab="",outside.points=F,pointcol=tfcols[i],main=expression(U < 0 ~ "vs" ~ beta > 0 ))
    for( i in rev(tf)){
        par(new=T)
        hdrcde::hdr.boxplot.2d(nonneutral.posteriors$qc[[i]]$betaDistrib.new.neg,nonneutral.posteriors$qc[[i]]$utility.new.2,prob=prob,shadecols=tfcols[i],xlim=c(0,1),ylim=c(0,1),xlab="% negatively polarized individual",ylab="% of positively polarized rumors",outside.points=F,pointcol=tfcols[i])
    }

    hdr.boxplot.2d(nonneutral.parameters$betaDistrib.new.pos,nonneutral.parameters$utility.new,prob=prob,shadecols=alpha("grey",.9),xlim=c(0,1),ylim=c(0,1)  ,xlab="",ylab="",outside.points=F,pointcol=tfcols[i])
    for( i in rev(tf)){
        par(new=T)
        hdrcde::hdr.boxplot.2d(nonneutral.posteriors$qc[[i]]$betaDistrib.new.pos,nonneutral.posteriors$qc[[i]]$utility.new,prob=prob,shadecols=tfcols[i],xlim=c(0,1),ylim=c(0,1),xlab="% positively polarized individual",ylab="% of negatively polarized rumors",outside.points=F,pointcol=tfcols[i])
    }


    par(mfrow=c(1,2))
    hdr.boxplot.2d(nonneutral.parameters$betaDistrib.new.neg,nonneutral.parameters$betaDistrib.new.pos,prob=prob,shadecols=alpha("grey",.9),xlim=c(0,1),ylim=c(0,1)  ,xlab="",ylab="",outside.points=F,pointcol=tfcols[i],main=expression(U < 0 ~ "vs" ~ beta > 0 ))
    for( i in rev(tf)){
        par(new=T)
        hdrcde::hdr.boxplot.2d(nonneutral.posteriors$qc[[i]]$betaDistrib.new.neg,nonneutral.posteriors$qc[[i]]$betaDistrib.new.pos,prob=prob,shadecols=tfcols[i],xlim=c(0,1),ylim=c(0,1),xlab="% negatively polarized individual",ylab="% of positively polarized rumors",outside.points=F,pointcol=tfcols[i])
    }

    hdr.boxplot.2d(nonneutral.parameters$betaDistrib.new.pos,nonneutral.parameters$utility.new,prob=prob,shadecols=alpha("grey",.9),xlim=c(0,1),ylim=c(0,1)  ,xlab="",ylab="",outside.points=F,pointcol=tfcols[i])
    for( i in rev(tf)){
        par(new=T)
        hdrcde::hdr.boxplot.2d(nonneutral.posteriors$qc[[i]]$betaDistrib.new.pos,nonneutral.posteriors$qc[[i]]$utility.new,prob=prob,shadecols=tfcols[i],xlim=c(0,1),ylim=c(0,1),xlab="% positively polarized individual",ylab="% of negatively polarized rumors",outside.points=F,pointcol=tfcols[i])
    }


    hdr.boxplot.2d(nonneutral.parameters[,5]+nonneutral.parameters[,6],nonneutral.parameters[,8]+nonneutral.parameters[,9],prob=prob,shadecols=alpha("grey",.9),xlim=c(0,1),ylim=c(0,1)  ,xlab="% - polarized individual",ylab="% - polarized individual",outside.points=F,pointcol=tfcols[i])
    for( i in rev(tf)){
        par(new=T)
        hdr.boxplot.2d(nonneutral.posteriors$qc[[i]][,5] + nonneutral.posteriors$qc[[i]][,6] ,nonneutral.posteriors$qc[[i]][,8] + nonneutral.posteriors$qc[[i]][,9],prob=prob,shadecols=tfcols[i],xlim=c(0,1),ylim=c(0,1),xlab="% polarized individual",ylab="% of polarized rumors",outside.points=F,pointcol=tfcols[i])
    }

    hdr.boxplot.2d(nonneutral.parameters$betaDistrib.new.neut,nonneutral.parameters$utility.new.pol,prob=prob,xlim=c(0,1),ylim=c(0,1) ,xlab="% neutral individual",ylab="% of polarized rumors",outside.points=F)
    for( i in (tf)){
        par(new=T)

        hdr.boxplot.2d(nonneutral.posteriors$qc[[i]]$betaDistrib.new.neut,nonneutral.posteriors$qc[[i]]$utility.new.pol,prob=prob,shadecol=tfcols[i],xlim=c(0,1),ylim=c(0,1),xlab="% neutral individual",ylab="% of polarized rumors",outside.points=F)
    }

}


printAllPosteriors <- function(){

    legnames=1:length(names(nonneutral.posteriors$qc$false))
    legnames=c("Number of active agents","Number of Rumors","timestep","Perceived % of tweets",expression(beta==-100),expression(beta==-10),expression(beta==0),expression(beta==10),expression(beta==100),expression(U==-1),expression(U==0),expression(U==1),expression(mu),"Initial number of Cascades","Total number of agents","decay time","backward time",expression(beta==0),expression(beta %in% group("{",list(-100,-10,10,100),"}")),expression(U==0),expression(U %in% group("{",list(-1,1),"}")),"ratio utilities (log(polarized/neutral))","ratio preferences (log(polarized/neutral))")
    names(legnames)=names(nonneutral.posteriors$qc$false)
    png("allposteriors_nonneutral.png",width=900,height=1200)
    par(mfrow=c(5,5))
    for(n in names(nonneutral.posteriors$qc$false))
        plot2dens(A=nonneutral.posteriors$qc$false[[n]],B=nonneutral.posteriors$qc$true[[n]],prior=nonneutral.parameters[[n]],xlab=legnames[n],main="",cols=cols)
    legend("topright",legend=c("prior","size","true"),fill=cols,cex=1)
    dev.off()

    par(mfrow=c(1,2))
    for(n in names(nonneutral.posteriors$qc$false)[22:23])
        plot2dens(A=nonneutral.posteriors$qc$false[[n]],B=nonneutral.posteriors$qc$true[[n]],prior=nonneutral.parameters[[n]],xlab=legnames[n],main="",cols=cols)


    par(mfrow=c(2,3))
    for(n in names(alberto.posteriors$qc$false))
        plot2dens(A=alberto.posteriors$qc$false[[n]],B=alberto.posteriors$qc$true[[n]],prior=alberto.parameters[[n]],xlab=n,main="",cols=cols)
    legend("topright",legend=c("prior","size","true"),fill=cols,cex=1)
    dev.new()

    par(mfrow=c(2,3))
    for(n in names(albertoc.posteriors$qc$false))
        plot2dens(A=albertoc.posteriors$qc$false[[n]],B=albertoc.posteriors$qc$true[[n]],prior=albertoc.parameters[[n]],xlab=n,main="",cols=cols)
    legend("topright",legend=c("prior","size","true"),fill=cols,cex=1)


    hdrcde::hdr.den(den=density(c(rnorm(100,5,1),rnorm(100,5,1))),prob = c(50,95),col=c("green","dark green"),bgcol="white",lwd=3,legend=T) 
    removeBadPrior <- function(){
        plot(density(alberto.parameters$mu * alberto.parameters$t_step * alberto.parameters$Nmin))
        plot(density(neutral.parameters$mu_c * neutral.parameters$repetition * neutral.parameters$IC))

        testscores=lapply(alberto.scores,lapply,function(u){u[alberto.parameters$mu < 0.0005] = 1000; return(u)})
        testalberto.posteriors=getAllposteriors(testscores,alberto.parameters,100)
        testscores.posteriors.scores.qc=lapply(testscores$qc,function(i)i[rank(i,ties="first")<100])

        tneutrscores=lapply(neutral.scores,lapply,function(u){u[neutral.parameters$mu < 0.0005] = 1000; return(u)})
        tneutrneutral.posteriors=getAllposteriors(tneutrscores,neutral.parameters,100)
        tneutrscores.posteriors.scores.qc=lapply(tneutrscores$qc,function(i)i[rank(i,ties="first")<100])
        lines(density(tneutrscores.posteriors.scores.qc$true))
        plot(density(testscores.posteriors.scores.qc$true))
        plot(density(neutral.parameters$mu_c * neutral.parameters$repetition * neutral.parameters$IC))

    }

    par(mfrow=c(2,3))
    for(n in names(post$topfiveAC$qc$false))
        plot2dens(A=post$topfiveAC$qc$false[[n]],B=post$topfiveAC$qc$true[[n]],prior=allparamaters$topfiveAC[[n]],xlab=n,main="",cols=cols)
    legend("topright",legend=c("prior","size","true"),fill=cols,cex=1)

    png("posteriors_neutral.png",width=480*1.4,height=480*1.6,pointsize=20)
    par(mfrow=c(2,2))
    par(mar=c(4,4,1,1))
    plot2dens(A=post[["neutral"]]$qc$false$mu,B=post[["neutral"]]$qc$true$mu,xlab="mu",main="",from=0,to=.002,cols=cols)
    legend("topright",legend=c("prior","false","true"),fill=cols,cex=1)
    plot2dens(A=post[["neutral"]]$qc$false$Nmin,B=post[["neutral"]]$qc$true$Nmin,prior=c(1000,10000),xlab="Nmin",main="",cols=cols)
    plot2dens(A=post[["neutral"]]$qc$false$t_step,B=post[["neutral"]]$qc$true$t_step,prior=c(10,300),xlab="t_step",main="",cols=cols)
    plot2dens(A=post[["neutral"]]$qc$false$tau,B=post[["neutral"]]$qc$true$tau,prior=c(1,100),xlab="tau",main="",cols=cols)
    dev.off()

}


