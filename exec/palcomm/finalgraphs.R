##Function that should centralize all script use to do write figure in Carrignon et al 2019
## but it still contain some "test"  that maybe be  usefull. Should be cleaned
library(devtools) 
load_all("../../")
library("hdrcde") #better to use hdrce last version

divide <- function(u)sapply(u,"/",u)

####GLOBAL VARIABLES
allru=c()
allru$size=tapply(allca$size,allca$rumor_id,sum)

trueru=c()
trueru$size=tapply(trueca$size,trueca$rumor_id,sum)

falseru=c()
falseru$size=tapply(falseca$size,falseca$rumor_id,sum)

mixedru=c()
mixedru$size=tapply(mixedca$size,mixedca$rumor_id,sum)

size=list()
size$true=trueru$size
size$false=falseru$size
tfcols=c("green","red")
names(tfcols)=c("true","false")
tf=c("true","false")
names(tf)=c("true","false")
listexp=c(
          "../testAllvsTFsinC/",
          "../testAllvsTF",
          "../testAllvsTFMine",
          "../testAllvsTFAlbertoCorrected/"
          #"../testneutralRumorsVF/"
          ) #this is use with getAllParameters but not use anymore as the folder have been read already and the results  save in `data_visualisation/`
names(listexp)=c("random","conformist","topfiveS","topfiveAC")#,"topfive3D")
########


GETALLDATA <- function(){


    ##first we get all score for all experiment
    alllims=lapply(listexp,checkIfNumberGood)
    allresults=lapply(names(listexp),function(u)getAllscores(listexp[[u]],lim=alllims[[u]],idscores = c("qc"),metrics=c("true","false"))) #can be long 
    ##Then all parameters correspondign to those score
    allparameters=lapply(names(listexp),function(u)getAllparameters(listexp[[u]],lim=alllims[[u]])) #can be long 
    ##Then we check the experiment with less simulation and will resample the result to have this number of everyone.
    names(allresults)=names(allparameters)=names(listexp)
    nmin=min(as.vector(sapply(allresults,function(u)lengths(u$qc)))) #to get the experiments with the smallest number of simulation.
    nmin=8770560
    #save(allresults,file="data_visualisation/all_getscore.bin")
    #save(allparameters,file="data_visualisation/all_parameters.bin")
    #load(file="data_visualisation/all_getscore.bin")
    for(exp in names(listexp)){
        subn=sample.int(length(allparameters[[exp]][[1]]),nmin)
        allresults[[exp]]=lapply(allresults[[exp]],lapply,"[",subn)
        allparameters[[exp]]=lapply(allparameters[[exp]],"[",subn)
    }

   # nmin=min(as.vector(sapply(allresults,function(u)lengths(u$qc)))) #to get the experiments with the smallest number of simulation.
   # allresults=lapply(allresults,lapply,lapply,function(u)u [sample(length(u),nmin)]) #this resampling is not good if we want to do it in phase with parameters in order to get the right posteriors
    #save(allresults,file="data_visualisation/normalized_all_getscore.bin")
    #save(allparameters,file="data_visualisation/normalized_all_parameters.bin")
    load(file="data_visualisation/normalized_all_getscore.bin")
    load(file="data_visualisation/normalized_all_parameters.bin")

    post=lapply(names(listexp),function(n)getAllposteriors(allresults[[n]],as.data.frame(allparameters[[n]]),1000))
    names(post)=names(listexp)
    save(post,file="data_visualisation/normalized_all_posteriors.bin")
    load(file="data_visualisation/normalized_all_posteriors.bin")
    rm(allresults,allparameters)
    #save(allparameters,file="data_visualisation/normalized_all_parameters.bin")
}

RERUNFROMPOSTERIORS <- function(post){

   rerun=list()
    for(exp in names(post)[4]){
        print(exp)
        cur=post[[exp]][[1]]
        alberto=F
        if(exp=="topfiveA")alberto=T
        if(exp=="topfiveAC")alberto=T 
        rerun[[exp]]=lapply(tf,function(l)reruns(cur[[l]],modelwrapper=randomCascades.list,data=size[[l]],scorefun=quantilediff,repet=0,samplepost=10000,type="sampling",alberto=alberto,log=F))

    }
save(rerun[[exp]],file="data_visualisation/back_exptopfiveAC")

}


HDRGRAPHS <- function(listexp){

    bins=10^(seq(0,5.5,.1))
    for( n in names(listexp)){
        print(n) ###in this part we use stored data and remove it right after it as the stored data are the full distribution of 10 000 simulation  + scores twice for each model which is pretty heavy ~ 15 - 24 Mo which led the memory to be full and R to crash
        load(paste0("data_visualisation/back_exp",n))  
        #exp=rerun[[n]]
        exp=lapply(exp,function(u)u[unlist(lapply(u,lapply,function(i)if(is.na(i$score)) F else T))]) #remove NA 
        allbined=lapply(exp,lapply,function(l)as.numeric(prop.table(table(cut(l[[1]]$distrib,breaks=bins,include.lowest = T)))))
        rm(exp)
        gc()
        tiff(paste0("hdr_check_",n,".tiff"),width=1200,height=600,pointsize=15)
        par(mfrow=c(1,2))
        ##Should Write this in a function and move to the package in ABC.R:
        na=lapply(tf,function(d)
                  {
                      allb=allbined[[d]]
                      par(mar=c(4,4,1,1))
                      dfallbin=do.call("rbind",allb)
                      dfallbin[dfallbin==0]=NA
                      dfallbin.log=log10(dfallbin)
                      hdr.boxplot(lapply((1:ncol(dfallbin))[-c(2,3,6)],function(i){
                                                 alog=dfallbin.log[,i]
                                                 alog[is.na(alog)]=-7
                                                 alog
}
                      ),space = .1,col=alpha(tfcols[d],c(.1,.3)),ylim=c(-7.2,0),prob = c(50,95),yaxt="n",outline = T,lambda=1,h=.1,pch=".",xlab="Number of RTs",ylab="Frequencies",from=0)
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T)))))[-c(2,3,6)],col=tfcols[d],pch=20)
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T)))))[-c(2,3,6)])
                      axis(1,at = (1:length(bins))[0:6*10+1],labels =  bins[0:6*10+1])
                      axis(2,at = -7:0,labels =  c(0,10^(-6:0)))
                      legend(1,-5.2,lty=c(NA,NA,1,NA),pch=c(NA,NA,NA,21),fill=c(alpha(tfcols[d],c(.1,.3)),NA,NA),border=c(1,1,NA,NA),legend=c("95% HDR","50% HDR","mode","data"),pt.bg=c(NA,NA,NA,tfcols[d]))
                  }

        )
        dev.off()
    }

    load("data_visualisation/back_exprandom")
    allbined=lapply(exp,lapply,function(l)as.numeric(prop.table(table(cut(l[[1]]$distrib,breaks=bins,include.lowest = T)))))
    truedist=as.numeric(prop.table(table(cut(size[["true"]],breaks=bins,include.lowest = T))))

    allbined=do.call("rbind",allbined$true)
    png("images/onebinDensity.png",width=600,height=600,pointsize=15)
    par(mar=c(4,4,1,1))
    hdr.den(den=density(allbined[, 31]),xlab="% of simulations",ylab="density",prob = c(50,95),col=alpha(tfcols["true"],c(.2,.5)),main="",bgcol=0,lwd=3,plot.line=F,legend=T) 
    #polygon(c(truedist[31]-.0002,truedist[31]+.0002,truedist[31]+.0002,truedist[31]-.0002),c(0.1,0.1,35,35),col="green",lwd=.5,border=1)
    points(truedist[31],.6,cex=1,pch=25,bg="green")
    text(truedist[31],.6,cex=.9,label="data",pos=4)
    #legend("right",pch=c(NA,NA,25),legend=c("50% HDR","95% HDR","data"),fill=c(alpha(tfcols["true"],c(.2,.5)),NA),pt.bg=c(NA,NA,"green"),border=c(1,1,NA))
    dev.off()


}

bayesFactorCalculation <- function(){

    load("data_visualisation/distribayesfactor.bin")
    load("data_visualisation/all_getscore.bin")
    modeldistrib=getModelDistribution(allresults)
    rm(allresults)

    total=length(modeldistrib$true$data)
    per=c(100,10,1)*500/total #if we want the top500/5000/50000
    cats=per*total
    names(cats)=paste(per*100,"%")
    ##to determin the epsilon 
    ##0.0001 * total
    ##0.00001 * total
    ##0.00002 * total
    ##0.000015 * total
    nparam=c('conformist'=3,'random'=2,'topfiveS'=4,'topfiveAC'=4)
    ##compute likelhood and AIC
    AIC=lapply(m,function(i)-2*log(i[,3]*cats[3]/nmin)+2*nparam) 
    m=lapply(modeldistrib,function(l)sapply(cats,function(c)table(l$model[rank(l$data,ties.method="min")<=c])))
    #save(m,file="data_visualisation/distribution_model_given_cats_withWrongAlberto.bin")
    #save(modeldistrib,file="data_visualisation/distribution_models_withWrongAlberto.bin")
    #load(file="data_visualisation/distribution_models_withWrongAlberto.bin")

    modeldistrib=lapply(modeldistrib,function(ml)list(data=ml$data[ml$model != "topfiveA"],model=ml$model[ml$model != "topfiveA"]))

    total=length(modeldistrib$true$data)
    cats=round(per*total)
    names(cats)=paste(per*100,"%")
    load(file="data_visualisation/distribution_models.bin")
    m=lapply(modeldistrib,function(l)sapply(cats,function(c)table(l$model[rank(l$data,ties.method="min")<=c])))
    ##save(m,file="data_visualisation/distribution_model_given_cats.bin")
    load(file="data_visualisation/distribution_model_given_cats.bin")
    ##save(modeldistrib,file="data_visualisation/distribution_models.bin")
    m=lapply(m,function(u)apply(u,2,function(i)round(i/sum(i),digit=3))) #get the ratio
    mBF=lapply(m,function(u)round(divide(u[,3]),digit=2))


    ###subsample model distribution given nonneutral.scores
    ##remove prior to far from nonneutral prior
    for(exp in names(listexp)){
        allresults[[exp]]=lapply(allresults[[exp]],lapply,"[",(allparameters[[exp]]$mu < 0.035 & allparameters[[exp]]$t_step > 50))
        allparameters[[exp]]=allparameters[[exp]][allparameters[[exp]]$mu < 0.035 & allparameters[[exp]]$t_step > 50,]
    }

    nmin=length(nonneutral.scores$qc$true)
    #reduce number of simulation to the same
    
    per=c(100,10,1)*500/nmin
    mfalse01=c()
    mfalse1=c()
    mfalse001=c()
    for(i in 1:50){
    suballresults=allresults
    suballparameters=allparameters
    for(exp in names(listexp)){
        subn=sample.int(length(allparameters[[exp]][[1]]),nmin)
        suballresults[[exp]]=lapply(allresults[[exp]],lapply,"[",subn)
        suballparameters[[exp]]=lapply(allparameters[[exp]],"[",subn)
    }

    undersample=getModelDistribution(suballresults)
    mdl="content"
    undersample$true$data=c(undersample$true$data,nonneutral.scores$qc$true)
    undersample$false$data=c(undersample$false$data,nonneutral.scores$qc$false)
    undersample$true$model=c(undersample$true$model,rep(mdl,lenundr))
    undersample$false$model=c(undersample$false$model,rep(mdl,lenundr))

    total=length(undersample$true$data)
    cats=round(per*total)
    names(cats)=paste(per*100,"%")

    m=lapply(undersample,function(l)sapply(cats,function(c)table(factor(l$model[rank(l$data,ties.method="min")<=c],levels=c("conformist","content","random","topfiveAC","topfiveS")))))
    m=lapply(m,function(u)apply(u,2,function(i)round(i/sum(i),digit=3))) #get the ratio
    mfalse01=rbind(mfalse01,m$false[,2])
    mfalse001=rbind(mfalse001,m$false[,3])
    mfalse1=rbind(mfalse1,m$false[,1])
    }
    mBF=lapply(m,function(u)round(divide(u[,2]),digit=2))

    colnames(mfalse01)=c("Conformist","Content","Unbiased","Top Alberto","Top Threshold")
    colnames(mfalse001)=c("Conformist","Content","Unbiased","Top Alberto","Top Threshold")
    pdf("images/top50content.pdf",pointsize=15,width=9,height=5)
    boxplot(mfalse01,ylab="% of the top 50",main="")
    dev.off()
    pdf("images/top5content.pdf",pointsize=15,width=9,height=5)
    boxplot(mfalse001,ylab="% of the top 5",main="")
    dev.off()
    #save(modeldistrib,file="data_visualisation/distribution_models_withcontent.bin")
    load(file="data_visualisation/distribution_models_withcontent.bin")
    for(mdl in names(listexp)){
        tmp=modeldistrib$true$data[modeldistrib$true$model == mdl]
        lenundr=nmin
        undr=sample.int(length(tmp),lenundr)
        undersample$true$data=c(undersample$true$data,tmp[undr])
        undersample$false$data=c(undersample$false$data,tmp[undr])
        undersample$true$model=c(undersample$true$model,rep(mdl,lenundr))
        undersample$false$model=c(undersample$false$model,rep(mdl,lenundr))
    } #nonneutral.score is computed in analysesABC.R
    for(mdl in names(listexp)){
        tmp=modeldistrib$true$data[modeldistrib$true$model == mdl]
        lenundr=nmin
        undr=sample.int(length(tmp),lenundr)
        undersample$true$data=c(undersample$true$data,tmp[undr])
        undersample$false$data=c(undersample$false$data,tmp[undr])
        undersample$true$model=c(undersample$true$model,rep(mdl,lenundr))
        undersample$false$model=c(undersample$false$model,rep(mdl,lenundr))
    } #nonneutral.score is computed in analysesABC.R
}

plotHDRnonneutral <- function(){
    bins=2^(seq(0,21,1))
        allbined=lapply(rerun[["content"]],lapply,function(l)as.numeric(prop.table(table(cut(l,breaks=bins,include.lowest = T)))))
       # png(paste0("hdr_check_6bins",n,".png"),width=1200,height=600,pointsize=15)
        par(mfrow=c(1,2))
        ##Should Write this in a function and move to the package in ABC.R:
        na=lapply(tf,function(d)
                  {
                      allb=allbined[[d]]
                      par(mar=c(4,4,1,1),cex=1.7)
                      dfallbin=do.call("rbind",allb)
                      dfallbin[dfallbin==0]=NA
                      dfallbin.log=log10(dfallbin)
                      ymin=-4.2
                      hdr.boxplot(lapply((1:ncol(dfallbin))[-c(2,3,6)],function(i){
                                                 alog=dfallbin.log[,i]
                                                 alog[is.na(alog)]=-100
                                                 alog
                      }
                      ),space = .1,col=alpha(tfcols[d],c(.1,.3)),ylim=c(ymin,0),prob = c(50,95),yaxt="n",outline = T,lambda=1,h=.1,pch=".",xlab="Number of RTs",ylab="Frequencies")
                      counts=lapply((1:ncol(dfallbin)),function(i){alog=dfallbin.log[,i];length(alog[alog == -100])})
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T))))),bg=tfcols[d],pch=21)
                      lbls=sapply(seq(0,length(bins),length.out = 5),function(i)as.expression(bquote(2^.(round(i)))))
                      axis(1,at = seq(1,length(bins),length.out = 5),labels = lbls)
                      axis(2,at = -3:0,labels =  c(10^(-3:0)))
                      #legend(1,ymin+2,lty=c(NA,NA,1,NA),pch=c(NA,NA,NA,21),fill=c(alpha(tfcols[d],c(.1,.3)),NA,NA),border=c(1,1,NA,NA),legend=c("95% HDR","50% HDR","mode","data"),pt.bg=c(NA,NA,NA,tfcols[d]))
                      legend("topright",lty=c(NA,NA,1,NA),pch=c(NA,NA,NA,21),fill=c(alpha(tfcols[d],c(.1,.3)),NA,NA),border=c(1,1,NA,NA),legend=c("95% HDR","50% HDR","mode","data"),pt.bg=c(NA,NA,NA,tfcols[d]),cex=.8)
                      #lines(0:length(bins),round(unlist($true)/10000,digits=2))
                      #text(,rep(-4.5,length(bins)),,cex=.5)
                      lines(1:1000,rep(ymin,1000),lty=2)
                      lines(1:1000,rep(ymin+1/3,1000),lty=2)
                      text(2,ymin-0.02,"0%",cex=.55,pos=2)
                      text(2,ymin-0.02+1/3+0.05,"100%",cex=.55,pos=2)
                      text(2,ymin-0.02+1/5,"% of sim. ended with 0 cascades of this size",cex=.55,pos=4)
                      lines(ymin+(round(unlist(counts)/length(allbined$false),digits=2))/3,col=tfcols[d])
                      list(zeros=counts,bins=allb)
                  }

        )
        dev.off()
    }


plotAllPosteriors <- fuinction(){
   cols=c(alpha("grey",.9),alpha("red",.8),alpha("green",.8))
    tiff("images/posteriors_random.tiff",width=480*1.4,height=480*1.6,pointsize=20)
    par(mfrow=c(2,2))
    par(mar=c(4,2,1,1),cex=1.05)
    plot2dens(A=post[["random"]]$qc$false$mu,B=post[["random"]]$qc$true$mu,prior=c(0,.3),xlab=expression(mu),main="",from=0,to=.003,cols=cols,yaxt="n")
    legend("topright",legend=c("prior","false","true"),fill=cols,cex=1)
    plot2dens(A=post[["random"]]$qc$false$Nmin,B=post[["random"]]$qc$true$Nmin,prior=c(1000,10000),xlab=expression("N"),main="",cols=cols,yaxt="n")
    plot2dens(A=post[["random"]]$qc$false$t_step,B=post[["random"]]$qc$true$t_step,prior=c(10,300),xlab=expression(t[step]),main="",cols=cols,yaxt="n")
    plot2dens(A=post[["random"]]$qc$false$tau,B=post[["random"]]$qc$true$tau,prior=c(1,100),xlab=expression(tau),main="",cols=cols,yaxt="n")
    dev.off()

    pdf("images/2_posteriors_random.pdf",width=8*1.4,height=8*1.6,pointsize=18)
    par(mfrow=c(2,2))
    par(mar=c(4,2,1,1),cex=1.05)
    plot2dens(A=post[["random"]]$qc$false$mu,B=post[["random"]]$qc$true$mu,prior=c(0,.3),xlab=expression(mu),main="",from=0,to=.003,cols=cols,yaxt="n")
    legend("topright",legend=c("prior","false","true"),fill=cols,cex=1)
    plot2dens(A=post[["random"]]$qc$false$Nmin,B=post[["random"]]$qc$true$Nmin,prior=c(1000,10000),xlab=expression("N"),main="",cols=cols,yaxt="n")
    plot2dens(A=post[["random"]]$qc$false$t_step,B=post[["random"]]$qc$true$t_step,prior=c(10,300),xlab=expression(t[step]),main="",cols=cols,yaxt="n")
    plot2dens(A=post[["random"]]$qc$false$tau,B=post[["random"]]$qc$true$tau,prior=c(1,100),xlab=expression(tau),main="",cols=cols,yaxt="n")
    dev.off()

    png("images/posteriors_conformist.png",width=480*1.4,height=480*1.8,pointsize=20)
    par(mfrow=c(3,2))
    par(mar=c(4,4,1,1))
    plot2dens(A=post[["conformist"]]$qc$false$mu,B=post[["conformist"]]$qc$true$mu,xlab="mu",main="",prior=c(0,.002),cols=cols)
    legend("topright",legend=c("prior","false","true"),fill=cols,cex=1)
    plot2dens(A=post[["conformist"]]$qc$false$Nmin,B=post[["conformist"]]$qc$true$Nmin,prior=c(1000,10000),xlab="Nmin",main="",cols=cols)
    plot2dens(A=post[["conformist"]]$qc$false$t_step,B=post[["conformist"]]$qc$true$t_step,prior=c(10,300),xlab="t_step",main="",cols=cols)
    plot2dens(A=post[["conformist"]]$qc$false$tau,B=post[["conformist"]]$qc$true$tau,prior=c(1,100),xlab="tau",main="",cols=cols)
    plot2dens(A=post[["conformist"]]$qc$false$beta,B=post[["conformist"]]$qc$true$beta,prior=c(-2,2),xlab="beta",main="",cols=cols)
    dev.off()


    png("images/posteriors_topfiveS.png",width=480*1.4,height=480*3,pointsize=20)
    par(mfrow=c(4,2))
    par(mar=c(4,4,1,1))
    plot2dens(A=post[["topfiveS"]]$qc$false$mu,B=post[["topfiveS"]]$qc$true$mu,xlab="mu",main="",prior=c(0,0.3),cols=cols)
    legend("topright",legend=c("prior","false","true"),fill=cols,cex=1)
    plot2dens(A=post[["topfiveS"]]$qc$false$Nmin,B=post[["topfiveS"]]$qc$true$Nmin,prior=c(1000,10000),xlab="Nmin",,main="",cols=cols)
    plot2dens(A=post[["topfiveS"]]$qc$false$t_step,B=post[["topfiveS"]]$qc$true$t_step,prior=c(10,300),xlab="t_step",main="",cols=cols)
    plot2dens(A=post[["topfiveS"]]$qc$false$tau,B=post[["topfiveS"]]$qc$true$tau,prior=c(1,100),xlab="tau",main="",cols=cols)
    plot2dens(A=post[["topfiveS"]]$qc$false$C,B=post[["topfiveS"]]$qc$true$C,prior=c(0,1),xlab="C",main="",cols=cols)
    plot2dens(A=post[["topfiveS"]]$qc$false$TF,B=post[["topfiveS"]]$qc$true$TF,prior=c(0,1000),xlab="tf",main="",cols=cols)
    plot2dens(A=log10( post$topfiveS$qc$false$TF/post$topfiveS$qc$false$Nmin ),col=cols,B=log10( post$topfiveS$qc$true$TF/post$topfiveS$qc$true$Nmin ),prior=log10(sample(1000:10000,100000,replace=T)/sample(100,100000,replace=T) ),from=-5,to=5,xlab="top  as % of pop == log(TF/NMIN)",main="")
    plot2dens(A=log10( post$topfiveS$qc$false$TF/post$topfiveS$qc$false$Nmin * post$topfiveS$qc$false$C ),col=cols,B=log10( post$topfiveS$qc$true$TF/post$topfiveS$qc$true$Nmin * post$topfiveS$qc$true$C),prior=log10(sample(1000:10000,100000,replace=T)/sample(100,100000,replace=T) * runif(100000) ),from=-5,to=5,xlab="% of people to copy == log10(TF/NMIN * C)",main="")
    dev.off()

    
    png("images/posteriors_topfiveA.png",width=480*1.4,height=480*3,pointsize=20)
    par(mfrow=c(4,2))
    par(mar=c(4,4,1,1))
    plot2dens(A=post[["topfiveA"]]$qc$false$mu,B=post[["topfiveA"]]$qc$true$mu,xlab="mu",main="",from=0,to=0.01,cols=cols)
    legend("topright",legend=c("prior","false","true"),fill=cols,cex=1)
    plot2dens(A=post[["topfiveA"]]$qc$false$Nmin,B=post[["topfiveA"]]$qc$true$Nmin,prior=c(1000,10000),xlab="Nmin",main="",cols=cols)
    plot2dens(A=post[["topfiveA"]]$qc$false$t_step,B=post[["topfiveA"]]$qc$true$t_step,prior=c(10,300),xlab="t_step",main="",cols=cols)
    plot2dens(A=post[["topfiveA"]]$qc$false$tau,B=post[["topfiveA"]]$qc$true$tau,prior=c(1,100),xlab="tau",main="",cols=cols)
    plot2dens(A=post[["topfiveA"]]$qc$false$C,B=post[["topfiveA"]]$qc$true$C,prior=c(0,1),xlab="C",main="",cols=cols)
    plot2dens(A=post[["topfiveA"]]$qc$false$TF,B=post[["topfiveA"]]$qc$true$TF,prior=c(0,1000),xlab="tf",main="",cols=cols)
    plot2dens(A=log10( post$topfiveA$qc$false$TF/post$topfiveA$qc$false$Nmin ),col=cols,B=log10( post$topfiveA$qc$true$TF/post$topfiveA$qc$true$Nmin ),prior=log10(sample(1000:10000,100000,replace=T)/sample(100,100000,replace=T) ),from=-4,to=4,xlab="top  as % of pop == log10(TF/NMIN)",main="")
    plot2dens(A=log10( post$topfiveA$qc$false$TF/post$topfiveA$qc$false$Nmin * post$topfiveA$qc$false$C ),col=cols,B=log10( post$topfiveA$qc$true$TF/post$topfiveA$qc$true$Nmin * post$topfiveA$qc$true$C),prior=log10(sample(1000:10000,100000,replace=T)/sample(100,100000,replace=T) * runif(100000) ),from=-4,to=4,xlab="% of people to copy == log10(TF/NMIN * C)",main="")
    dev.off()
    dev.new()

    png("images/posteriors_topfiveAC.png",width=480*1.4,height=480*3,pointsize=20)
    par(mfrow=c(4,2))
    par(mar=c(4,4,1,1))
    plot2dens(A=post[["topfiveAC"]]$qc$false$mu,B=post[["topfiveAC"]]$qc$true$mu,xlab="mu",main="",from=0,to=0.01,cols=cols)
    legend("topright",legend=c("prior","false","true"),fill=cols,cex=1)
    plot2dens(A=post[["topfiveAC"]]$qc$false$Nmin,B=post[["topfiveAC"]]$qc$true$Nmin,prior=c(1000,10000),xlab="Nmin",main="",cols=cols)
    plot2dens(A=post[["topfiveAC"]]$qc$false$t_step,B=post[["topfiveAC"]]$qc$true$t_step,prior=c(10,300),xlab="t_step",main="",cols=cols)
    plot2dens(A=post[["topfiveAC"]]$qc$false$tau,B=post[["topfiveAC"]]$qc$true$tau,prior=c(1,100),xlab="tau",main="",cols=cols)
    plot2dens(A=post[["topfiveAC"]]$qc$false$C,B=post[["topfiveAC"]]$qc$true$C,prior=c(0,1),xlab="C",main="",cols=cols)
    plot2dens(A=post[["topfiveAC"]]$qc$false$TF,B=post[["topfiveAC"]]$qc$true$TF,prior=c(0,1000),xlab="tf",main="",cols=cols)
    plot2dens(A=log10( post$topfiveAC$qc$false$TF/post$topfiveAC$qc$false$Nmin ),col=cols,B=log10( post$topfiveAC$qc$true$TF/post$topfiveAC$qc$true$Nmin ),prior=log10(sample(1000:10000,100000,replace=T)/sample(100,100000,replace=T) ),from=-4,to=4,xlab="top  as % of pop == log10(TF/NMIN)",main="")
    plot2dens(A=log10( post$topfiveAC$qc$false$TF/post$topfiveAC$qc$false$Nmin * post$topfiveAC$qc$false$C ),col=cols,B=log10( post$topfiveAC$qc$true$TF/post$topfiveAC$qc$true$Nmin * post$topfiveAC$qc$true$C),prior=log10(sample(1000:10000,100000,replace=T)/sample(100,100000,replace=T) * runif(100000) ),from=-4,to=4,xlab="% of people to copy == log10(TF/NMIN * C)",main="")
    dev.off()
    dev.new()
}


getAllHdrAndDoTablle <- function(post){
    myhdr <- function(a)unlist(hdr(a,lambda=.000001,prob=c(95))[c("mode","hdr")])

    for(m in names(post)){
        try({allhdr=lapply(post[[m]],lapply,lapply,myhdr)
        alltable=lapply(allhdr,do.call,what=rbind)
        for(u in colnames(alltable$qc))print(file=paste0("tables/interval_model_",m,"_param_",u,".tex"),xtable(do.call(rbind,alltable$qc[,u]),caption=sanitize(paste("mode and 95 % interval for the paramater",u, "of the model",m)),digits=NULL,auto=T,label=paste("tab","hdr",m,u,sep=":")))
        })
        print(m)
    }

}

allCDF <- function(){
    for( n in names(listexp)[1:4]){
        print(n) ###in this part we use stored data and remove it right after it as the stored data are the full distribution of 10 000 simulation  + scores twice for each model which is pretty heavy ~ 15 - 24 Mo which led the memory to be full and R to crash
        load(paste0("data_visualisation/back_exp",n))  
        #par(mfrow=c(1,2))
        ##Should Write this in a function and move to the package in ABC.R:
        png(paste0("cdf_check_",n,".png"),width=1200,height=600,pointsize=15)
        par(mfrow=c(1,2))
        na=lapply(tf,function(d)
                  {
                      plotCCFD(size[[d]],col=tfcols[d],main=n)
                      lapply(exp[[d]][sample.int(1000)],function(u)pointsCCFD(u[[1]]$distrib,col=alpha(tfcols[d],.01)))
                      pointsCCFD(size[[d]],col=1,pch=1)
                  }

        )
        dev.off()
        rm(exp)

    }

}

randomHDR  <- function(){
    bins=10^(seq(0,6,1))
    for( n in names(listexp)[1:2]){
        print(n) ###in this part we use stored data and remove it right after it as the stored data are the full distribution of 10 000 simulation  + scores twice for each model which is pretty heavy ~ 15 - 24 Mo which led the memory to be full and R to crash
        load(paste0("data_visualisation/back_exp",n))  
        #exp=rerun[[n]]
        exp=lapply(exp,function(u)u[unlist(lapply(u,lapply,function(i)if(is.na(i$score)) F else T))]) #remove NA 
        allbined=lapply(exp,lapply,function(l)as.numeric(prop.table(table(cut(l[[1]]$distrib,breaks=bins,include.lowest = T)))))
        rm(exp)
        gc()
        png(paste0("hdr_check_6bins",n,".png"),width=1200,height=600,pointsize=15)
        par(mfrow=c(1,2))
        ##Should Write this in a function and move to the package in ABC.R:
        na=lapply(tf,function(d)
                  {
                      allb=allbined[[d]]
                      par(mar=c(4,4,1,1))
                      dfallbin=do.call("rbind",allb)
                      dfallbin[dfallbin==0]=NA
                      dfallbin.log=log10(dfallbin)
                      hdr.boxplot(lapply((1:ncol(dfallbin)),function(i){
                                         alog=dfallbin.log[,i]
                                         alog[is.na(alog)]=-7
                                         alog
                      }
                      ),space = .1,col=alpha(tfcols[d],c(.1,.3)),ylim=c(-7.2,0),prob = c(50,95),yaxt="n",outline = T,lambda=1,h=.1,pch=".",xlab="Number of RTs",ylab="Frequencies",from=0)
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T))))),col=tfcols[d],pch=20)
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T))))))
                      axis(1,at = (1:length(bins)),labels =  bins)
                      axis(2,at = -7:0,labels =  c(0,10^(-6:0)))
                      legend(1,-5.2,lty=c(NA,NA,1,NA),pch=c(NA,NA,NA,21),fill=c(alpha(tfcols[d],c(.1,.3)),NA,NA),border=c(1,1,NA,NA),legend=c("95% HDR","50% HDR","mode","data"),pt.bg=c(NA,NA,NA,tfcols[d]))
                  }

        )
        dev.off()
    }


    bins=2^(seq(0,6,1))
    for( n in names(listexp)[1:2]){
        print(n) ###in this part we use stored data and remove it right after it as the stored data are the full distribution of 10 000 simulation  + scores twice for each model which is pretty heavy ~ 15 - 24 Mo which led the memory to be full and R to crash
        load(paste0("data_visualisation/back_exp",n))  
        #exp=rerun[[n]]
        exp=lapply(exp,function(u)u[unlist(lapply(u,lapply,function(i)if(is.na(i$score)) F else T))]) #remove NA 
        allbined=lapply(exp,lapply,function(l)as.numeric(prop.table(table(cut(l[[1]]$distrib,breaks=bins,include.lowest = T)))))
        rm(exp)
        gc()
        png(paste0("hdr_check_6bins",n,".png"),width=1200,height=600,pointsize=15)
        par(mfrow=c(1,2))
        ##Should Write this in a function and move to the package in ABC.R:
        na=lapply(tf,function(d)
                  {
                      allb=allbined[[d]]
                      par(mar=c(4,4,1,1))
                      dfallbin=do.call("rbind",allb)
                      dfallbin[dfallbin==0]=NA
                      dfallbin.log=log10(dfallbin)
                      hdr.boxplot(lapply((1:ncol(dfallbin)),function(i){
                                         alog=dfallbin.log[,i]
                                         alog[is.na(alog)]=-7
                                         alog
                      }
                      ),space = .1,col=alpha(tfcols[d],c(.1,.3)),ylim=c(-7.2,0),prob = c(50,95),yaxt="n",outline = T,lambda=1,h=.1,pch=".",xlab="Number of RTs",ylab="Frequencies",from=0)
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T))))),col=tfcols[d],pch=20)
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T))))))
                      axis(1,at = (1:length(bins)),labels =  bins[0:7])
                      axis(2,at = -7:0,labels =  c(0,10^(-6:0)))
                      legend(1,-5.2,lty=c(NA,NA,1,NA),pch=c(NA,NA,NA,21),fill=c(alpha(tfcols[d],c(.1,.3)),NA,NA),border=c(1,1,NA,NA),legend=c("95% HDR","50% HDR","mode","data"),pt.bg=c(NA,NA,NA,tfcols[d]))
                  }

        )
        dev.off()
    }

    bins=2^(seq(0,17,1))
    for( n in names(listexp)){
        print(n) ###in this part we use stored data and remove it right after it as the stored data are the full distribution of 10 000 simulation  + scores twice for each model which is pretty heavy ~ 15 - 24 Mo which led the memory to be full and R to crash
        load(paste0("data_visualisation/back_exp",n))  
        exp=lapply(exp,function(u)u[unlist(lapply(u,lapply,function(i)if(is.na(i$score)) F else T))]) #remove NA 
        allbined=lapply(exp,lapply,function(l)as.numeric(prop.table(table(cut(l[[1]]$distrib,breaks=bins,include.lowest = T)))))
        rm(exp)
        gc()
        png(paste0("hdr_check_power2bins",n,".png"),width=1200,height=600,pointsize=15)
        par(mfrow=c(1,2))
        ##Should Write this in a function and move to the package in ABC.R:
        na=lapply(tf,function(d)
                  {
                      allb=allbined[[d]]
                      par(mar=c(4,4,1,1))
                      dfallbin=do.call("rbind",allb)
                      dfallbin[dfallbin==0]=NA
                      dfallbin.log=log10(dfallbin)
                      hdr.boxplot(lapply((1:ncol(dfallbin)),function(i){
                                         alog=dfallbin.log[,i]
                                         alog[is.na(alog)]=-7
                                         alog
                      }
                      ),space = .1,col=alpha(tfcols[d],c(.1,.3)),ylim=c(-7.2,0),prob = c(50,95),yaxt="n",outline = T,lambda=1,h=.1,pch=".",xlab="Number of RTs",ylab="Frequencies",from=0)
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T))))),col=tfcols[d],pch=20)
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T))))))
                      axis(1,at = (1:length(bins)),labels =  bins)
                      axis(2,at = -7:0,labels =  c(0,10^(-6:0)))
                      legend(1,-5.2,lty=c(NA,NA,1,NA),pch=c(NA,NA,NA,21),fill=c(alpha(tfcols[d],c(.1,.3)),NA,NA),border=c(1,1,NA,NA),legend=c("95% HDR","50% HDR","mode","data"),pt.bg=c(NA,NA,NA,tfcols[d]))
                  }

        )
        dev.off()
    }


    bins=10^(seq(0,5.5,.1))
    for( n in names(listexp)[5]){
        print(n) ###in this part we use stored data and remove it right after it as the stored data are the full distribution of 10 000 simulation  + scores twice for each model which is pretty heavy ~ 15 - 24 Mo which led the memory to be full and R to crash
        exp=rerun[[n]]
        exp=lapply(exp,function(u)u[unlist(lapply(u,lapply,function(i)if(is.na(i$score)) F else T))]) #remove NA 
        allbined=lapply(exp,lapply,function(l)as.numeric(prop.table(table(cut(l[[1]]$distrib,breaks=bins,include.lowest = T)))))
        rm(exp)
        gc()
        png(paste0("hdr_check_cut0_",n,".png"),width=1200,height=600,pointsize=15)
        par(mfrow=c(1,2))
        ##Should Write this in a function and move to the package in ABC.R:
        na=lapply(tf,function(d)
                  {
                      allb=allbined[[d]]
                      par(mar=c(4,4,1,1))
                      dfallbin=do.call("rbind",allb)
                      dfallbin[dfallbin==0]=NA
                      dfallbin.log=log10(dfallbin)
                      hdr.boxplot(lapply((1:ncol(dfallbin))[-c(2,3,6)],function(i){
                                         alog=dfallbin.log[,i]
                                         alog[is.na(alog)]=-8
                                         alog
                      }
                      ),space = .1,col=alpha(tfcols[d],c(.1,.3)),ylim=c(-7.2,0),prob = c(50,95),yaxt="n",outline = T,lambda=1,h=.1,pch=".",xlab="Number of RTs",ylab="Frequencies",from=0)
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T)))))[-c(2,3,6)],col=tfcols[d],pch=20)
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T)))))[-c(2,3,6)])
                      axis(1,at = (1:length(bins))[0:6*10+1],labels =  bins[0:6*10+1])
                      axis(2,at = -6:0,labels =  c(10^(-6:0)))
                      legend(1,-5.2,lty=c(NA,NA,1,NA),pch=c(NA,NA,NA,21),fill=c(alpha(tfcols[d],c(.1,.3)),NA,NA),border=c(1,1,NA,NA),legend=c("95% HDR","50% HDR","mode","data"),pt.bg=c(NA,NA,NA,tfcols[d]))
                  }

        )
        dev.off()
    }

    bins=10^(seq(0,5.5,.1))
    for( n in names(listexp)){
        print(n) ###in this part we use stored data and remove it right after it as the stored data are the full distribution of 10 000 simulation  + scores twice for each model which is pretty heavy ~ 15 - 24 Mo which led the memory to be full and R to crash
        load(paste0("data_visualisation/back_exp",n))  
        exp=lapply(exp,function(u)u[unlist(lapply(u,lapply,function(i)if(is.na(i$score)) F else T))]) #remove NA 
        allbined=lapply(exp,lapply,function(l)as.numeric(prop.table(table(cut(l[[1]]$distrib,breaks=bins,include.lowest = T)))))
        rm(exp)
        gc()
        png(paste0("hdr_check_no0_",n,".png"),width=1200,height=600,pointsize=15)
        par(mfrow=c(1,2))
        ##Should Write this in a function and move to the package in ABC.R:
        ymin=-5.2
        na=lapply(tf,function(d)
                  {
                      allb=allbined[[d]]
                      par(mar=c(4,4,1,1))
                      dfallbin=do.call("rbind",allb)
                      dfallbin[dfallbin==0]=NA
                      dfallbin.log=log10(dfallbin)
                      hdr.boxplot(lapply((1:ncol(dfallbin))[-c(2,3,6)],function(i){
                                         alog=dfallbin.log[,i]
                                         alog[is.na(alog)]=-100
                                         alog
                      }
                      ),space = .1,col=alpha(tfcols[d],c(.1,.3)),ylim=c(ymin,0),prob = c(50,95),yaxt="n",outline = T,lambda=1,h=.1,pch=".",xlab="Number of RTs",ylab="Frequencies")
                      counts=lapply((1:ncol(dfallbin))[-c(2,3,6)],function(i){alog=dfallbin.log[,i];length(alog[alog == -100])})
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T)))))[-c(2,3,6)],col=tfcols[d],pch=20)
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T)))))[-c(2,3,6)])
                      axis(1,at = (1:length(bins))[0:6*10+1],labels =  bins[0:6*10+1])
                      axis(2,at = -6:0,labels =  c(10^(-6:0)))
                      legend(1,ymin+2,lty=c(NA,NA,1,NA),pch=c(NA,NA,NA,21),fill=c(alpha(tfcols[d],c(.1,.3)),NA,NA),border=c(1,1,NA,NA),legend=c("95% HDR","50% HDR","mode","data"),pt.bg=c(NA,NA,NA,tfcols[d]))
                      #text((1:length(bins)),rep(-4.5,length(bins)),round(unlist(na$true)/10000,digits=2),cex=.5)
                  }

        )
        dev.off()
    }

    bins=10^(seq(0,5.5,.1))
    for( n in names(listexp)[1]){
        print(n) ###in this part we use stored data and remove it right after it as the stored data are the full distribution of 10 000 simulation  + scores twice for each model which is pretty heavy ~ 15 - 24 Mo which led the memory to be full and R to crash
        load(paste0("data_visualisation/back_exp",n))  
        exp=lapply(exp,function(u)u[unlist(lapply(u,lapply,function(i)if(is.na(i$score)) F else T))]) #remove NA 
        allbined=lapply(exp,lapply,function(l)as.numeric(prop.table(table(cut(l[[1]]$distrib,breaks=bins,include.lowest = T)))))
        rm(exp)
        gc()
        tiff(paste0("hdr_check_count0_",n,".tiff"),width=1200,height=600,pointsize=15)
        par(mfrow=c(1,2))
        ##Should Write this in a function and move to the package in ABC.R:
        na=lapply(tf,function(d)
                  {
                      allb=allbined[[d]]
                      par(mar=c(4,4,1,1))
                      dfallbin=do.call("rbind",allb)
                      dfallbin[dfallbin==0]=NA
                      dfallbin.log=log10(dfallbin)
                      ymin=-5.2
                      hdr.boxplot(lapply((1:ncol(dfallbin))[-c(2,3,6)],function(i){
                                         alog=dfallbin.log[,i]
                                         alog[is.na(alog)]=-100
                                         alog
                      }
                      ),space = .1,col=alpha(tfcols[d],c(.1,.3)),ylim=c(ymin,0),prob = c(50,95),yaxt="n",outline = T,lambda=1,h=.1,pch=".",xlab="Number of RTs",ylab="Frequencies")
                      counts=lapply((1:ncol(dfallbin))[-c(2,3,6)],function(i){alog=dfallbin.log[,i];length(alog[alog == -100])})
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T)))))[-c(2,3,6)],col=tfcols[d],pch=20)
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T)))))[-c(2,3,6)])
                      axis(1,at = (1:length(bins))[0:6*10+1],labels =  bins[0:6*10+1])
                      axis(2,at = -6:0,labels =  c(10^(-6:0)))
                      legend(1,ymin+2,lty=c(NA,NA,1,NA),pch=c(NA,NA,NA,21),fill=c(alpha(tfcols[d],c(.1,.3)),NA,NA),border=c(1,1,NA,NA),legend=c("95% HDR","50% HDR","mode","data"),pt.bg=c(NA,NA,NA,tfcols[d]))
                      #lines(0:length(bins),round(unlist($true)/10000,digits=2))
                      #text(,rep(-4.5,length(bins)),,cex=.5)
                      lines(ymin+(round(unlist(counts)/10000,digits=2))/3,col=tfcols[d])
                      abline(h=ymin,lty=2)
                      abline(h=ymin+1/3,lty=2)
                      text(1,ymin-0.02,"0%",cex=.6)
                      text(1,ymin-0.02+1/3+0.05,"100%",cex=.6)
                      text(1,ymin-0.02+1/5,"% of sim. ended with 0 cascades of this size",cex=.6,pos=4)
                  }

        )
        dev.off()
    }

    bins=10^(seq(0,5.5,.1))
    for( n in names(listexp)){
        print(n) ###in this part we use stored data and remove it right after it as the stored data are the full distribution of 10 000 simulation  + scores twice for each model which is pretty heavy ~ 15 - 24 Mo which led the memory to be full and R to crash
        load(paste0("data_visualisation/back_exp",n))  
        exp=lapply(exp,function(u)u[unlist(lapply(u,lapply,function(i)if(is.na(i$score)) F else T))]) #remove NA 
        allbined=lapply(exp,lapply,function(l)as.numeric(prop.table(table(cut(l[[1]]$distrib,breaks=bins,include.lowest = T)))))
        rm(exp)
        gc()
        png(paste0("hdr_check_nolog_",n,".png"),width=1200,height=600,pointsize=15)
        ymin=0
        par(mfrow=c(1,2))
        ##Should Write this in a function and move to the package in ABC.R:
        na=lapply(tf,function(d)
                  {
                      allb=allbined[[d]]
                      par(mar=c(4,4,1,1))
                      dfallbin=do.call("rbind",allb)
                      hdr.boxplot(lapply((1:ncol(dfallbin))[-c(2,3,6)],function(u)dfallbin[,u]),space = .1,col=alpha(tfcols[d],c(.1,.3)),prob = c(50,95),outline = T,h=.1,pch=".",xlab="Number of RTs",ylab="Frequencies",ylim=c(0,max(dfallbin)))
                      points(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T))))[-c(2,3,6)],col=tfcols[d],pch=20)
                      points(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T))))[-c(2,3,6)])
                      axis(1,at = (1:length(bins))[0:6*10+1],labels =  bins[0:6*10+1])
                  }

        )
        dev.off()
    }
}

    THEBINUSEINTHEPAPER <- function(){
    bins=2^(seq(0,21,1))
    for( n in names(listexp)[3:4]){
        print(n) ###in this part we use stored data and remove it right after it as the stored data are the full distribution of 10 000 simulation  + scores twice for each model which is pretty heavy ~ 15 - 24 Mo which led the memory to be full and R to crash
        load(paste0("data_visualisation/back_exp",n))  
        exp=lapply(exp,function(u)u[unlist(lapply(u,lapply,function(i)if(is.na(i$score)) F else T))]) #remove NA 
        if(n=="content")
            allbined=lapply(rerun[["content"]],lapply,function(l)as.numeric(prop.table(table(cut(l,breaks=bins,include.lowest = T))))) #i didnt get rerun from the rerun function
        else
            allbined=lapply(exp,lapply,function(l)as.numeric(prop.table(table(cut(l[[1]]$distrib,breaks=bins,include.lowest = T)))))
        rm(exp)
        gc()
        #tiff(paste0("hdr_check_count0_",n,".tiff"),width=1200,height=600,pointsize=15)
        pdf(paste0("hdr_check_count0_",n,".pdf"),width=18,height=9,pointsize=15)
        par(mfrow=c(1,2))
        ##Should Write this in a function and move to the package in ABC.R:
        na=lapply(tf,function(d)
                  {
                      allb=allbined[[d]]
                      par(mar=c(4,4,1,1),cex=1.7)
                      dfallbin=do.call("rbind",allb)
                      dfallbin[dfallbin==0]=NA
                      dfallbin.log=log10(dfallbin)
                      ymin=-4.2
                      hdr.boxplot(lapply((1:ncol(dfallbin))[-c(2,3,6)],function(i){
                                                 alog=dfallbin.log[,i]
                                                 alog[is.na(alog)]=-100
                                                 alog
                      }
                      ),space = .1,col=alpha(tfcols[d],c(.1,.3)),ylim=c(ymin,0),prob = c(50,95),yaxt="n",outline = T,lambda=1,h=.1,pch=".",xlab="Number of RTs",ylab="Frequencies")
                      counts=lapply((1:ncol(dfallbin)),function(i){alog=dfallbin.log[,i];length(alog[alog == -100])})
                      points(log10(as.numeric(prop.table(table(cut(size[[d]],breaks=bins,include.lowest = T))))),bg=tfcols[d],pch=21)
                      lbls=sapply(seq(0,length(bins),length.out = 5),function(i)as.expression(bquote(2^.(round(i)))))
                      axis(1,at = seq(1,length(bins),length.out = 5),labels = lbls)
                      axis(2,at = -3:0,labels =  c(10^(-3:0)))
                      #legend(1,ymin+2,lty=c(NA,NA,1,NA),pch=c(NA,NA,NA,21),fill=c(alpha(tfcols[d],c(.1,.3)),NA,NA),border=c(1,1,NA,NA),legend=c("95% HDR","50% HDR","mode","data"),pt.bg=c(NA,NA,NA,tfcols[d]))
                      legend("topright",lty=c(NA,NA,1,NA),pch=c(NA,NA,NA,21),fill=c(alpha(tfcols[d],c(.1,.3)),NA,NA),border=c(1,1,NA,NA),legend=c("95% HDR","50% HDR","mode","data"),pt.bg=c(NA,NA,NA,tfcols[d]),cex=.8)
                      #lines(0:length(bins),round(unlist($true)/10000,digits=2))
                      #text(,rep(-4.5,length(bins)),,cex=.5)
                      lines(1:1000,rep(ymin,1000),lty=2)
                      lines(1:1000,rep(ymin+1/3,1000),lty=2)
                      text(2,ymin-0.02,"0%",cex=.55,pos=2)
                      text(2,ymin-0.02+1/3+0.05,"100%",cex=.55,pos=2)
                      text(2,ymin-0.02+1/5,"% of sim. ended with 0 cascades of this size",cex=.55,pos=4)
                      lines(ymin+(round(unlist(counts)/10000,digits=2))/3,col=tfcols[d])
                      list(zeros=counts,bins=allb)
                  }

        )
        dev.off()
    }
    }




    #########
    graphRumorsNew <- function(){
        alldate=as.POSIXct(allca$date)
        allca$hours=round(alldate,"hours") #trunc by hours
        allca$days=round(alldate,"days")     #trunc by days
        allca$month=round(as.POSIXct(paste0(strftime(alldate,format="%Y-%m"),"-01"), tz="GMT"),"days")     #trunc by month
        allca$years=round(as.POSIXct(paste0(strftime(alldate,format="%Y"),"-01-01"), tz="GMT"),"days")     #trunc by year 
        par(mfrow=c(1,2),mar=c(5,5,1,1))
        rateintro=sort(tapply(allca$month-min(allca$month),allca$rumor_id,min)) #for each rumor we assign the date of the first cascade that's spread it.  alldayz=seq(min(rateintro),max(rateintro),3600*24) #generate all possible month, will be used to avoid skipping the day whith not new rumors in order to to take into account the day where number of rumor = 0
        rateintro=as.POSIXct(rateintro,origin=min(allca$month))
        plot(table(rateintro))
        plot(unique(rateintro),as.numeric(table(rateintro)),type="h",lwd=5,xlab="",ylab="Number of new rumors")


    }


    graphScoreRobustness <- function(){
        lotofnmin=sapply(seq(1,20,2),function(s)c(s,summary(apply(post$random$qc$true[sample.int(length(post$random$qc$true[[1]]),1000),],1,function(i){print(s);i["Nmin"] = i["Nmin"] * s;quantilediff(randomCascades.list(unlist(i))$size,trueru$size)}))))
        lotof=sapply(seq(1,10,1.5),function(s)c(s,summary(apply(post$random$qc$true[sample.int(length(post$random$qc$true[[1]]),500),],1,function(i){print(s);i["t_step"] = i["t_step"] * s;quantilediff(randomCascades.list(unlist(i))$size,trueru$size)}))))


    }


    dataorubt <- function(){

        tiff("images/1a_cascades_size.tiff",width=600,height=600,pointsize=20)
        par(mar=c(4,4,1,1))
        plotCCFD(allca$size,main="",axes=F,xlim=c(1,100000),ylim=c(0.0001,100),xlab="cascades size",type="n")
    #pointsCCFD(mixedca$size,col="orange")
    pointsCCFD(falseca$size,col="red")  
    pointsCCFD(trueca$size ,col="green")
        axis(1,at=10^(0:5),label=sapply(0:5,function(i)as.expression(bquote(10^.(i)))))
        axis(2,at=c(0.0001,0.01,1,100),label=c(expression(10^-4),expression(10^-2),1,100))
 legend("bottomleft",pch=c(20,20),col=c("red","green"),legend=c("false","true") )
        box()
        dev.off()

        pdf("images/1a_cascades_size.pdf",width=9,height=9,,pointsize=18)
        par(mar=c(4,4,1,1))
        plotCCFD(allca$size,main="",axes=F,xlim=c(1,100000),ylim=c(0.0001,100),xlab="cascades size",type="n")
    #pointsCCFD(mixedca$size,col="orange")
    pointsCCFD(falseca$size,col="red")  
    pointsCCFD(trueca$size ,col="green")
        axis(1,at=10^(0:5),label=sapply(0:5,function(i)as.expression(bquote(10^.(i)))))
        axis(2,at=c(0.0001,0.01,1,100),label=c(expression(10^-4),expression(10^-2),1,100))
 legend("bottomleft",pch=c(20,20),col=c("red","green"),legend=c("false","true") )
        box()
        dev.off()

        
        tiff("images/1b_aggregate_size.tiff",width=600,height=600,pointsize=20)
        par(mar=c(4,4,1,1))
        plotCCFD(allru$size,main="",axes=F,xlim=c(1,200000),ylim=c(0.03,100),xlab="cascades size",type="n")
    #pointsCCFD(mixedru$size,col="orange")
    pointsCCFD(falseru$size,col="red")  
    pointsCCFD(trueru$size ,col="green")
        axis(1,at=10^(0:5),label=sapply(0:5,function(i)as.expression(bquote(10^.(i)))))
        axis(2,at=c(0.1,1,10,100),label=c(expression(10^-1),1,10,100))
 legend("bottomleft",pch=c(20,20),col=c("red","green"),legend=c("false","true") )
        box()
        dev.off()

        pdf("images/1b_aggregate_size.pdf",width=9,height=9,pointsize=18)
        par(mar=c(4,4,1,1))
        plotCCFD(allru$size,main="",axes=F,xlim=c(1,200000),ylim=c(0.03,100),xlab="cascades size",type="n")
    #pointsCCFD(mixedru$size,col="orange")
    pointsCCFD(falseru$size,col="red")  
    pointsCCFD(trueru$size ,col="green")
        axis(1,at=10^(0:5),label=sapply(0:5,function(i)as.expression(bquote(10^.(i)))))
        axis(2,at=c(0.1,1,10,100),label=c(expression(10^-1),1,10,100))
 legend("bottomleft",pch=c(20,20),col=c("red","green"),legend=c("false","true") )
        box()
        dev.off()
    }


    ###SIMPL3D run and plot
    #woule2=(getRumors(cascades3D(N=400000,R=100000,IC=1000,Nmax=2000,betadistrib=rep(0,400000),utility=rep(0,100000),metrics="size",time=200,mu_c=.00002,memory=F,captl=5000,dtime=0,stime=0)))

    #plotCCFD(allru$size,main="",axes=F,xlim=c(1,200000),ylim=c(0.03,100),xlab="cascades size",type="n")
    ##pointsCCFD(mixedru$size,col="orange")
    #pointsCCFD(falseru$size,col="red")  
    #pointsCCFD(trueru$size ,col="green")
    #axis(1,at=10^(0:5),label=sapply(0:5,function(i)as.expression(bquote(10^.(i)))))
    #axis(2,at=c(0.1,1,10,100),label=c(expression(10^-1),1,10,100))
    #legend("bottomleft",pch=c(20,20,20),col=c("red","green",1),legend=c("false","true","content bias (U=0,beta=0)") )
    #box()
    #pointsCCFD(woule2)
    #dev.off()

    # for(i in 1:10){
    #     load(paste0("reruncontent_",i));
    #     rerun$content$false=c(rerun$content$false,reruncontent$false)
    #     rerun$content$true=c(rerun$content$true,reruncontent$true)
    # }
    getModelDistribution <- function(allresults){
        modeldistrib=list()
        modeldistrib$true$data=lapply(allresults,function(u)u$qc$true)
        modeldistrib$false$data=lapply(allresults,function(u)u$qc$false)
        modeldistrib$true$model=rep(names(modeldistrib$true$data),lengths(modeldistrib$true$data))
        modeldistrib$false$model=rep(names(modeldistrib$false$data),lengths(modeldistrib$false$data))
        modeldistrib$true$data=unlist(lapply(allresults,function(u)u$qc$true),use.names = F)
        modeldistrib$false$data=unlist(lapply(allresults,function(u)u$qc$false),use.names = F)
        return(modeldistrib)
    }
