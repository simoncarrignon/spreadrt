#' From a main folder get all scores 
#' In theory, the folder `mainfold` should contain the results of an abc as done by ./abcdir/abc.R
#' which, in theory again, should look like a set of subfolder `mainfoldresults/mainfoldresultsX` where X depends on the number of simulation and the way they are divided in abc.R
#'@param mainfold mainfolder with allscores 
#'@param idscores the id of the score
#'@param metrics the ids of the metrics use for the score
#'@param log if true, show the name of the folder scanned
#'@return a list of list of list
#'
#'@export getAllscores 
getAllscores <- function(mainfold,idscores=c("kc","kokc","qc","kcrev","kokcrev"),metrics=c("depth","breadth","size"),log=T,lim=NULL){
    #With all the folder we can now get all the scores
    if(is.character(lim))
        allfolders=lim
    else
        allfolders=getlistfoldsubfold(mainfold,lim=lim)
    names(idscores)=idscores
    names(metrics)=metrics
    #That take a while. use of Rmpi could be good.
    allscores=lapply(idscores,function(iscore)
                     {
                         tis=lapply(allfolders,function(j)
                                    {
                                        if(log)print(j)
                                        load(file.path(j,"scores.bin"));
                                        if(length(metrics)<2)
                                            sapply(scores,function(tmps)
                                                   tryCatch(unlist(tmps[[iscore]][metrics]),error=function(e){res=rep(NA,length(metrics));names(res)=metrics;res})
                                                   )
                                        else
                                            as.data.frame(
                                                          t(
                                                            sapply(scores,function(tmps)
                                                                   tryCatch(unlist(tmps[[iscore]][metrics]),error=function(e){res=rep(NA,length(metrics));names(res)=metrics;res})
                                                                   )
                                                            )
                                                          )
                                    })
                         if(length(metrics)<2){
                             res=c()
                             res[[metrics]]=unlist(tis)
                             tis=res
                         }
                         else
                             tis=do.call("rbind",tis)
                         tis
                     }
    )
    names(allscores)=idscores
    return(allscores)
}

##This function is a hack when dimension of metrics and idscores == 1, that SHOULD be included in previous function as a test case
#'
#' @export getAllscoresB 
getAllscoresB <- function(mainfold,idscores=c("kc","kokc","qc","kcrev","kokcrev"),metrics=c("depth","breadth","size"),log=T,lim=NULL){
    #With all the folder we can now get all the scores
    if(is.character(lim))
        allfolders=lim
    else
        allfolders=getlistfoldsubfold(mainfold,lim=lim)
    names(idscores)=idscores
    names(metrics)=metrics
    #That take a while. use of Rmpi could be good.

    tis=sapply(allfolders,function(j){
               if(log)print(j)
               load(file.path(j,"scores.bin"));
               sapply(scores,function(at)as.numeric(at["qc"][[1]][[1]]))
    })
    res=c()
    res[[idscores]]=c()
    res[[idscores]][[metrics]]=tis
    return(res)
}

#' From a main folder get all parameters
#' In theory, the folder `mainfold` should contain the results of an abc as done by ./abcdir/abc.R
#' which, in theory again, should look like a set of subfolder `mainfoldresults/mainfoldresultsX` where X depends on the number of simulation and the way they are divided in abc.R
#'@param mainfold mainfolder with allscores 
#'@param log if true, show the name of the folder scanned
#'@param lim if lim is a vector of int then take the folder that have index in lim, if lim is a character, directly take lim as the list of folders
#'@param check if check is true then when check if all parameters led to good simulation (ie number of parameters -== number of score)
#' note that with check use of checkIfNumberGood should be useless
#'@return a dataframe with all parameters
#'
#' @export getAllparameters 
getAllparameters <- function(mainfold,log=F,lim=NULL,check=T){
    if(is.character(lim))
        allfolders=lim
    else
        allfolders=getlistfoldsubfold(mainfold,lim=lim)
    allparameters=lapply(allfolders,function(i)
                         {
                             if(log)print(i)
                             load(file.path(i,"parameters.bin"))
                             if(check){ ##check if all simulation went well
                                 load(file.path(i,"scores.bin"))
                                 o=length(scores)
                                 if(o != length(parameters[[1]])){ #if length of parameter is different than number of score we have to rely on the id stored to make the correspondance score/parameters:
                                    o=sapply(scores,function(u)tryCatch(u$id,error=function(e){NA})) #get ids of simulations that were stored w/ score
                                    if(log)print(paste("num ids:",length(o)))
                                    parameters=lapply(parameters,"[",o)
                                 }
                             }
                             return(as.data.frame(parameters))
                         })
    allparameters.dataframe=do.call("rbind.data.frame",allparameters) 
    return(allparameters.dataframe)
}


#' @export getlistfoldsubfold 
getlistfoldsubfold <- function(mainfold,lim){
    allfolders= list.dirs(mainfold,recursive=F)
    if(!is.null(lim))allfolders=paste0(file.path(mainfold,basename(mainfold)),lim)
    return(allfolders)
}


#' Return posterior given a priori distribution, score and a number of particle to accept
#'@param allscores a list of list of score
#'@param allparameters.dataframe a dataframe with all parameters
#'@param n an integer
#'@return a list of list, with for each score and each metrics a dataframe with subset of the parameter for with score of the simulation is in top `n`
#'
#' @export getAllposteriors 
getAllposteriors <- function(allscores,allparameters.dataframe,n=500)
    lapply(allscores,lapply,function(u)allparameters.dataframe[rank(u,ties="first")<n,])

#' Return the best set of parameter given a score an 
#'@param allscores a list of list of score
#'@param allparameters.dataframe a dataframe with all parameters
#'@param s id of the function use to compute the score
#'@param m metric use to compute the score
#'@return the list of parameter that gave the best score
#'
#' @export getbest 
getbest <- function(allscores,allparameters.dataframe,s,m)
    unlist(allparameters.dataframe[which(rank(allscores[[s]][[m]],ties="first") == 1),])

#' Wrapper for Cascades3D that automatically convert a list of arguments 
#'@param listofparameters a list of parameters normally given to cascades3D
#'@return same as cascades3D
#'
#' @export cascades3D.list 
cascades3D.list <- function(listofparameters,metrics=c("size","depth","breadth")){
    if(length(listofparameters)==9)
        cascades3D(
                   log=F,
                   N=listofparameters["N"],
                   stime=listofparameters["stime"],
                   dtime=listofparameters["dtime"],
                   IC=listofparameters["IC"],
                   Nmax=listofparameters["Nmax"],
                   mu_c=listofparameters["mu_c"],
                   R=listofparameters["R"],
                   time=listofparameters["repetition"],
                   captl=listofparameters["captl"],
                   utility = rep(0,listofparameters["R"]),
                   metrics=metrics,
                   betadistrib = rep(0,listofparameters["N"])
                   )
    else
    {
        prior_beta=c(-100,-10,0,10,100)
        prior_utility=c(-1,0,1)
        util=grep("utility.*",names(listofparameters))
        betad=grep("betaDistrib.*",names(listofparameters))
        cascades3D(
                   log=F,
                   metrics=metrics,
                   N=listofparameters["N"],
                   stime=listofparameters["stime"],
                   dtime=listofparameters["dtime"],
                   IC=listofparameters["IC"],
                   Nmax=listofparameters["Nmax"],
                   R=listofparameters["R"],
                   time=listofparameters["repetition"],
                   captl=listofparameters["captl"],
                   utility = generateDistribeFromFrequencies(prior_utility,listofparameters["R"],listofparameters[util]),
                   betadistrib = generateDistribeFromFrequencies(prior_beta,listofparameters["N"],listofparameters[betad])
                   )
    }
}

#' Wrapper for random cascade 
#'@param listofparameters a list of parameters normally given to cascades3D
#'@return same as randomCascades
#'
#' @export randomCascades.list 
randomCascades.list <- function(listofparameters,...){
    z=list(...)
    alberto=NULL
    if(!is.null(z$alberto))alberto=z$alberto

    if(length(listofparameters) == 6){ #topfive
        print(paste0("random cascade topfive alberto=",alberto))
    return(randomCascades(
                   Nmin=listofparameters["Nmin"],
                   mu=listofparameters["mu"],
                   t_step=listofparameters["t_step"],
                   tau=listofparameters["tau"],
                   conformity=F,
                   topfive=T,
                   C=listofparameters["C"],
                   TF=listofparameters["TF"],
                   ...
                   ))
    }
    if(length(listofparameters) == 5){#conformity
        print(paste0("random cascade conformity"))

    return(randomCascades(
                   Nmin=listofparameters["Nmin"],
                   mu=listofparameters["mu"],
                   t_step=listofparameters["t_step"],
                   tau=listofparameters["tau"],
                   conformity=T,
                   beta=listofparameters["beta"],
                   ...
                   ))
    }
    if(length(listofparameters) == 4){#random 

        print(paste0("random cascade neutral"))
    return(randomCascades(
                   Nmin=listofparameters["Nmin"],
                   mu=listofparameters["mu"],
                   t_step=listofparameters["t_step"],
                   tau=listofparameters["tau"],
                   topfive=F,
                   conformity=F,
                   ...
                   ))
    }
}




#' go within a folder with result and extract the good simulation
#'@param ind integer of the indice (ie id) of the simulation 
#'@param expefold folder within which we have to look
#'@param full boolean saying is the full simulation need to be retrieves
#'@return a list with the score and the list of parameters
#'
#' @export getSimuGivenIndices 
getSimuGivenIndices <- function(ind,expefold,full=F){
    i=1
    totest=0:(length(list.files("../testRumorSize/"))-1)
    for(t in totest){
        load(file.path(expefold,paste0(basename(expefold),t),"scores.bin"))
        for(s  in 1:length(scores)){
            if(ind == i){
                if(full){
                    load(file.path(expefold,paste0(basename(expefold),t),"fullsim",paste0("res_sim_",s,".bin")))
                    return(rud)
                }
                else{
                    load(file.path(expefold,paste0(basename(expefold),t),"parameters.bin"))
                    return(list(scores=scores[s],param=sapply(parameters,"[[",s)))
                }

            }
            i=i+1
        }
    }
}

#'Function that helps to run new test to check if score correspond do parameters.
#' allow to be sure that parameters and scores read from the folders make sense.
#' @param numberparam number of parameter that will be tested
#' @param rep number of time the parameter set will be repeted
#'
#' @export fromTestScoresVSParam 
fromTestScoresVSParam <- function(theparam,thescores,numberparam=20,repet=10,modelwrapper=cascades3D.list,data=allru$size,scorefun=quantilediff,rumor=F){
    stopifnot(length(thescores)== nrow(theparam))
    allsample=lapply(sample(length(thescores),numberparam),function(j)
                     {
                         nrep=replicate(repet,
                                        {
                                            simu=modelwrapper(unlist(theparam[j,]))
                                            if(rumor)rum=getRumors(simu)
                                            else rum=simu$size
                                            score=scorefun(rum,data)
                                            print(score)
                                            res=c()
                                            res$simu=simu
                                            res$score=score
                                            res$id=j
                                            res
                                        } )
                         return(nrep)
                     }
    )
    allsample.df=do.call("rbind.data.frame",lapply(allsample,function(at)as.data.frame(t(at[2:3,]))))
    allsample.df$score=unlist(allsample.df$score)
    allsample.df$id=unlist(allsample.df$id)
    allsample.df$simscore=thescores[allsample.df$id]
    allsample.df.mean=tapply(allsample.df$score,allsample.df$id,mean)
    allsample.df$imean=allsample.df.mean[as.character(allsample.df$id)]
    plot(allsample.df$score ~ allsample.df$imean,ylim=range(allsample.df$score,allsample.df$simscore))
    points(allsample.df$simscore ~ allsample.df$imean,col="red")
}

#' Same before but just return list of the simulation rerun 
#' @param numberparam number of parameter that will be tested
#' @param rep number of time the parameter set will be repeted
#'
#' @export reruns 
reruns <- function(theposterior,thescores=NULL,samplepost=100,repet=10,modelwrapper=cascades3D.list,data=allru$size,scorefun=quantilediff,rumor=F,type="all",log=F,...){
    ##To set the color scale I use the trick of calcualting all the simulation first and then attribue the color. Ungly so far.
    print(samplepost)
    #print(paste0("alberto=",alberto))
    if(type == "all" | type == "sampling"){
        if(type == "sampling"){
            resample=sample.int(nrow(theposterior),samplepost,replace=T)
            repet=1
        }
        else
            resample=sample.int(nrow(theposterior),samplepost)
        resampledpost=theposterior[resample,]
        if(!is.null(thescores)){
            resamplescore=thescores[resample]
            orderpost=resampledpost[order(resamplescore,decreasing=T),]
        }
        orderpost=resampledpost
    }
    else if(type == "median")
        orderpost=as.data.frame(lapply(theposterior,median))
    else if(type == "best")
        orderpost=as.data.frame(theposterior[which.max(thescores),])
    replpost=lapply(1:nrow(orderpost),function(r)
                   {
                       p=orderpost[r,]
                       return(lapply(1:repet,function(o)
                                      {
                                          if(log)print(o)
                                          sims=tryCatch(modelwrapper(unlist(p),...),error=function(e){print(e);NULL})
                                          if(is.null(sims))return(list(score=NA,distrib=NA))
                                          if(rumor) rums=getRumors(sims)
                                          else rums=sims$size
                                          sc=scorefun(rums,data)
                                          return(list(score=sc,distrib=rums))
                                      }
                       )
                       )
                   })
    return(replpost)
}


#' Same before but print distribution function that helps to run test to check if score correspond do parameters.
#' and then
#' @param numberparam number of parameter that will be tested
#' @param rep number of time the parameter set will be repeted
#' @param clrs if "score" the curves are colored givent the score function used 
#' @param ... parameters transmitted to plot
#'
#' @export plotPosteriorsCheck 
plotPosteriorsCheck <- function(theposterior,thescores,samplepost=100,repet=10,modelwrapper=cascades3D.list,data=allru$size,scorefun=quantilediff,rumor=F,type="all",clrs="score",alberto=F,...){

    reruns(theposterior,thescores,samplepost=100,repet=10,modelwrapper=modelwrapper,data=allru$size,scorefun=quantilediff,rumor=F,type="all",clrs="score",alberto=alberto)

    plotCCFD(data,...)
    if(clrs == "score"){
        print(replpost)
        allscores=unlist(lapply(replpost,lapply,"[[","score"))
        #rrs=range(allscores)
        #scale=(rrs[2]-rrs[1])/100
        #spacescore=seq(rrs[1],rrs[2],scale)
        #cols=topo.colors(length(spacescore))
        cols=topo.colors(length(unique(allscores)))
        names(cols)=sort(unique(allscores))
        na=lapply(replpost,lapply,function(s)pointsCCFD(s$distrib,col=alpha(cols[as.character(s$score)],.6)))

    }
    else
        lapply(replpost,lapply,function(s)pointsCCFD(s$distrib,col=clrs))
    return(replpost)
}

  
#' Take a folder name and return the list of id of the folders that have exact same number of parameters and scores
#' @param folder the name of a folder as the one used in getlistfoldsubfold
#'
#' @export checkIfNumberGood 
checkIfNumberGood <- function(folder,inv=T){
    alllength=lapply(getlistfoldsubfold(folder,lim=NULL),function(i){tryCatch({load(file.path(i,"parameters.bin"));load(file.path(i,"scores.bin"));c(i,length(parameters[[1]]),length(scores))},error=function(e)c(i,0,-1))})
    if(inv)return(unlist(lapply(1:length(alllength),function(i){s=alllength[[i]];if(s[2]==s[3])s[1]})))
    else return(unlist(lapply(1:length(alllength),function(i){s=alllength[[i]];if(s[2]!=s[3])s[1]})))
    return(rem)
}

