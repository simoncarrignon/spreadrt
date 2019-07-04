source("../spread/R/cascades3D.R")
source("../spread/R/tools.R")
source("../spread/R/ABC_analyse.R")

args=commandArgs(trailingOnly = TRUE) #pass number of slave by comand line, should be #node - 1 as one node is used by the master 

ns=args[1]#first argument is the number of slave
mainfold=args[2] #second argument = mainfolder to store the results
if(is.na(mainfold) | mainfold=="") mainfold="simulations"

print(paste0("Abc will be stored in mainfolder:",mainfold))
dir.create(mainfold)

print(getRumors)

library(Rmpi)
mpi.spawn.Rslaves(nslaves=ns)

mpi.bcast.cmd( id <- mpi.comm.rank() )
mpi.bcast.cmd( ns <- mpi.comm.size() )
mpi.bcast.cmd( host <- mpi.get.processor.name() )
mpi.bcast.cmd(source("../spread/R/cascades3D.R"))
mpi.bcast.cmd(source("../spread/R/tools.R"))
mpi.bcast.cmd(source("../spread/R/ABC_analyse.R"))
mpi.bcast.cmd(mainfold)
#(N=200,R=5,betadistrib=rep(1,200),utility=seq(0,1,length.out=5),repetition=2,tl=500,IC=50,summary=T,Nmax=NULL,lambda_c=0,lambda_r=0,dtime=0)

### Define Prior Ranges ###
prior_N=c(1000,60000)
prior_R=c(100,3000)
#prior_beta=c(-100,-10,0,10,100)
#prior_utility=c(-1,0,1)
prior_repetition=c(10,300)
prior_captl=c(1,600)
prior_mu_c=c(0,.3)
prior_IC=c(100,2500)
prior_Nmax=c(1000,10000)
prior_dtime=0:150
prior_stime=0:150

### Define Constants and Settings ###
nsim=768 #number of simulations
nsubfold=1
full=F #Boolean to save or not the full simulations

load("observations.bin")

fi=0
for(ns in 1:nsubfold){

    fold=file.path(mainfold,paste0(mainfold,fi))
    while(file.exists(fold)){
        fi=fi+1
        fold=file.path(mainfold,paste0(mainfold,fi))
    }
    dir.create(fold)
    if(full)dir.create(file.path(fold,"fullsim"))

    ### Create Parameter Space ###
    parameters=list(
                    ### Define Prior Ranges ###,
                    N=sample(prior_N[1]:prior_N[2],nsim,replace=T),
                    R=sample(prior_R[1]:prior_R[2],nsim,replace=T),
                    #betaDistrib=generalPartition(nsim,length(prior_beta)),
                    #betaDistrib_k=runif(nsim,prior_betaDistrib_k[1],prior_betaDistrib_k[2]),
                    #betaDistrib_t=runif(nsim,prior_betaDistrib_t[1],prior_betaDistrib_t[2]),
                    #utility=runif(nsim),
                    #utility_b=runif(nsim,prior_utility_b[1],prior_utility_b[2]),
                    #utility_a=runif(nsim,prior_utility_a[1],prior_utility_a[2]),
                    #utility=generalPartition(nsim,length(prior_utility)), #this will contain the percantage for each class of utility 
                    #utility=generalPartition(nsim,length(prior_utility)), #this will contain the percantage for each class of utility 
                    repetition=sample(prior_repetition[1]:prior_repetition[2],nsim,replace=T),
                    captl=runif(nsim,prior_captl[1],prior_captl[2]),
                    mu_c=runif(nsim,prior_mu_c[1],prior_mu_c[2]),
                    IC=sample(prior_IC[1]:prior_IC[2],nsim,replace=T),
                    Nmax=runif(nsim,prior_Nmax[1],prior_Nmax[2]),
                    dtime=sample(prior_dtime,nsim,replace=T),
                    stime=sample(prior_stime,nsim,replace=T)
                    )


    scores=mpi.applyLB(1:nsim,function(i,parameters,obs,fold,full){
                         start.time <- Sys.time()
                         param=lapply(parameters,"[[",i)
                         print(paste("BEFORE:",paste0(param,collapse=","),sep=","))
                         if(unlist(param)["mu_c"] < .1){
                             rud=cascades3D.list(unlist(param),metrics="size")

                             end.time <- Sys.time()
                             time.taken <- end.time - start.time
                             print(time.taken)
                             print(paste("mu was=",unlist(param)["mu_c"]))
                             print(paste("length rumors",length(getRumors(rud))))
                             print(paste("AFTER",paste0(param,collapse=","),as.numeric(difftime(start.time,end.time,units="secs")),sep=","))
                             rumsize=getRumors(rud)
                             if(full){
                                 tow=list(rd=rud,t=time.taken)
                                 #SBATCH --qos=debug
                                 save(tow,file=file.path(fold,"fullsim",paste0("res_sim_",i,".bin")))
                             }
                             qc=c()
                             print("calculate qc")
                             qc$size=quantilediff(rumsize , obs$allrsize)
                             qc$true=quantilediff(rumsize , obs$trueru$size)
                             qc$false=quantilediff(rumsize , obs$falseru$size)
                             print("done")
                             #qc$breadth=quantilediff(rud$breadth , obs$allbreadth)
                             #qc$depth=quantilediff(rud$depth , obs$alldepth+1)
                             kokc=c()
                             #kokc$size=kokldiff(rumsize , obs$allrsize)
                             #kokc$breadth=kokldiff(rud$breadth , obs$allbreadth)
                             #kokc$depth=kokldiff(rud$depth , obs$alldepth+1)
                             kokcrev=c()
                             #kokcrev$size=kokldiff(rumsize , obs$allrsize)
                             #kokcrev$breadth=kokldiff(obs$allbreadth , rud$breadth )
                             #kokcrev$depth=kokldiff(obs$alldepth+1,rud$depth)
                             kc=c()
                             #kc$breadth=kldiff(rud$breadth , obs$allbreadth)
                             #kc$depth=kldiff(rud$depth , obs$alldepth+1)
                             kcrev=c()
                             #kcrev$size=kldiff(rumsize , obs$allrsize)
                             #kcrev$breadth=kldiff(obs$allbreadth , rud$breadth )
                             #kcrev$depth=kldiff(obs$alldepth+1,rud$depth)
                             #print(list(kc,qc=qc,kokc=kokc,kcrev=kcrev,kokcrev=kokcrev))
                             print(list(kc,qc=qc,kokc=kokc,kcrev=kcrev,kokcrev=kokcrev))
                         }
                         else{
                             qc=c()
                             print("calculate qc")
                             qc$size=1000
                             qc$true=1000
                             qc$false=1000
                             print("done")
                             #qc$breadth=quantilediff(rud$breadth , obs$allbreadth)
                             #qc$depth=quantilediff(rud$depth , obs$alldepth+1)
                             kokc=c()
                             #kokc$size=kokldiff(rumsize , obs$allrsize)
                             #kokc$breadth=kokldiff(rud$breadth , obs$allbreadth)
                             #kokc$depth=kokldiff(rud$depth , obs$alldepth+1)
                             kokcrev=c()
                             #kokcrev$size=kokldiff(rumsize , obs$allrsize)
                             #kokcrev$breadth=kokldiff(obs$allbreadth , rud$breadth )
                             #kokcrev$depth=kokldiff(obs$alldepth+1,rud$depth)
                             kc=c()
                             #kc$breadth=kldiff(rud$breadth , obs$allbreadth)
                             #kc$depth=kldiff(rud$depth , obs$alldepth+1)
                             kcrev=c()
                             #kcrev$size=kldiff(rumsize , obs$allrsize)
                             #kcrev$breadth=kldiff(obs$allbreadth , rud$breadth )
                             #kcrev$depth=kldiff(obs$alldepth+1,rud$depth)
                             #print(list(kc,qc=qc,kokc=kokc,kcrev=kcrev,kokcrev=kokcrev))
                             print(list(kc,qc=qc,kokc=kokc,kcrev=kcrev,kokcrev=kokcrev))
                         }

                         #print(sapply(parameters,function(p)p[i]))
                         return(list(kc=kc,qc=qc,kokc=kokc,kcrev=kcrev,kokcrev=kokcrev,id=i))

                    },obs=observation,parameters=parameters,fold=fold,full=full)


    save(parameters,file=file.path(fold,"parameters.bin"))
    save(scores,file=file.path(fold,"scores.bin"))
}

mpi.close.Rslaves(dellog = T)
mpi.quit()

