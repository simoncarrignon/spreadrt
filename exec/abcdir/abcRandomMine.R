source("../spread/R/cascades3D.R")
source("../spread/R/tools.R")
source("../spread/R/ABC_analyse.R")

args=commandArgs(trailingOnly = TRUE) #pass number of slave by comand line, should be #node - 1 as one node is used by the master 

ns=args[1]#first argument is the number of slave
mainfold=args[2] #second argument = mainfolder to store the results
if(is.na(mainfold) | mainfold=="") mainfold="simulations"

print(paste0("Abc will be stored in mainfolder:",mainfold))
dir.create(mainfold)


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
prior_Nmin=c(1000,10000)
prior_mu=c(0,0.3)
#prior_Nmax=c(500,5000)
prior_t_step=c(10,300)
#prior_betaDistrib=c(1,5)
prior_C=c(0,1)
prior_TF=c(0,1000)
prior_tau=c(1,100)
#prior_beta=c(-2,2)










### Define Constants and Settings ###
nsim=7680 #number of simulations
nsubfold=50
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
    ## first ABC with conformity
    #parameters=list(
    #                ### Define Prior Ranges ###,
    #                Nmin=sample(prior_Nmin[1]:prior_Nmin[2],nsim,replace=T),
    #                Nmax=sample(prior_Nmax[1]:prior_Nmax[2],nsim,replace=T),
    #                mu=runif(nsim,prior_mu[1],prior_mu[2]),
    #                t_step=sample(prior_t_step[1]:prior_t_step[2],nsim,replace=T),
    #                tau=sample(prior_tau[1]:prior_tau[2],nsim,replace=T)
    #                )


    ### Create Parameter Space ###
    parameters=list(
                    ### Define Prior Ranges ###,
                    Nmin=sample(prior_Nmin[1]:prior_Nmin[2],nsim,replace=T),
                    mu=runif(nsim,prior_mu[1],prior_mu[2]),
                    t_step=sample(prior_t_step[1]:prior_t_step[2],nsim,replace=T),
                    tau=sample(prior_tau[1]:prior_tau[2],nsim,replace=T),
                    C=runif(nsim,prior_C[1],prior_C[2]),
                    TF=sample(prior_TF[1]:prior_TF[2],nsim,replace=T)
                    #beta=runif(nsim,prior_beta[1],prior_beta[2])
                    )












    scores=mpi.parLapply(1:nsim,function(i,parameters,obs,fold,full){
                         print(sapply(parameters,function(p)p[i]))
                         start.time <- Sys.time()
                         rumsize=randomCascades(
                                        Nmin=parameters$Nmin[i],
                                        mu=parameters$mu[i],
                                        conformity=F,
                                        t_step=parameters$t_step[i],
                                        tau=parameters$tau[i],
                                        C=parameters$C[i],
                                        TF=parameters$TF[i],
                                        alberto=F,
                                        topfive=T,
                                        beta=0
                                        )$size
                         end.time <- Sys.time()
                         time.taken <- end.time - start.time
                         print(time.taken)
                         if(full){
                             rud=list(scores=sc,rd=rumsize,t=time.taken)
                             save(rud,file=file.path(fold,"fullsim",paste0("res_sim_",i,".bin")))
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
                         print(list(kc=kc,qc=qc,kokc=kokc,kcrev=kcrev,kokcrev=kokcrev))
                         return(list(kc=kc,qc=qc,kokc=kokc,kcrev=kcrev,kokcrev=kokcrev,id=i))

                    },obs=observation,parameters=parameters,fold=fold,full=full)


    save(parameters,file=file.path(fold,"parameters.bin"))
    save(scores,file=file.path(fold,"scores.bin"))
}

mpi.close.Rslaves(dellog = T)
mpi.quit()

