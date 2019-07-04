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
prior_N=c(500,10000)
prior_R=c(10,1000)
prior_beta=c(-100,-10,0,10,100)
prior_utility=c(-1,0,1)
prior_repetition=c(5,80)
prior_captl=c(100,1000)
prior_lambda_c=c(0.00001,.1)
prior_IC=c(500,1000)
prior_Nmax=c(0.001,.1)
prior_dtime=c(1:5,-1)
prior_stime=c(1:5,-1)

### Define Constants and Settings ###
nsim=7680 #number of simulations
nsubfold=5
full=FALSE #Boolean to save or not the full simulations

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
                    betaDistrib=generalPartition(nsim,length(prior_beta)),
                    #betaDistrib_k=runif(nsim,prior_betaDistrib_k[1],prior_betaDistrib_k[2]),
                    #betaDistrib_t=runif(nsim,prior_betaDistrib_t[1],prior_betaDistrib_t[2]),
                    #utility=runif(nsim),
                    #utility_b=runif(nsim,prior_utility_b[1],prior_utility_b[2]),
                    #utility_a=runif(nsim,prior_utility_a[1],prior_utility_a[2]),
                    utility=generalPartition(nsim,length(prior_utility)), #this will contain the percantage for each class of utility 
                    repetition=sample(prior_repetition[1]:prior_repetition[2],nsim,replace=T),
                    captl=runif(nsim,prior_captl[1],prior_captl[2]),
                    lambda_c=runif(nsim,prior_lambda_c[1],prior_lambda_c[2]),
                    IC=sample(prior_IC[1]:prior_IC[2],nsim,replace=T),
                    Nmax=runif(nsim,prior_Nmax[1],prior_Nmax[2]),
                    dtime=sample(prior_dtime,nsim,replace=T),
                    stime=sample(prior_stime,nsim,replace=T)
                    )


    scores=mpi.applyLB(1:nsim,function(i,parameters,obs,fold,full){
                         print(sapply(parameters,function(p)p[i]))
                         start.time <- Sys.time()
                         rud=cascades3D(
                                        log=F,
                                        N=parameters$N[i],
                                        R=parameters$R[i],
                                        #betadistrib=rgamma(parameters$N[i],parameters$betaDistrib_k[i],parameters$betaDistrib_t[i]),
                                        betadistrib={
                                            b=rep(c(-100,-10,0,10,100),parameters$betaDistrib[i,]*parameters$N[i])
                                            if(sum(table(b))<parameters$N[i]){##In cas of missing value we ad to ad more, should be made a function, cannot think rn
                                                miss=parameters$N[i]-sum(table(b))
                                                b=c(b,sample(c(-100,-10,0,10,100),miss,repl=T))
                                                b
                                            }
                                            else b
                                        }
                                        ,
                                        #utility=rbeta(parameters$R[i],parameters$utility_a[i],parameters$utility_b[i]),
                                        utility={
                                            u=rep(c(-1,0,1),parameters$utility[i,]*parameters$R[i])
                                            if(sum(table(u))<parameters$R[i]){ ##In cas of missing value we ad to ad more
                                                miss=parameters$R[i]-sum(table(u))
                                                u=c(u,sample(c(-1,0,1),miss,repl=T))
                                                u
                                            }
                                            else u
                                        },
                                        time=parameters$repetition[i],
                                        captl=parameters$captl[i],
                                        lambda_c=parameters$lambda_c[i],
                                        Nmax=parameters$Nmax[i],
                                        dtime=parameters$dtime[i],
                                        stime=parameters$stime[i]
                                        )
                         end.time <- Sys.time()
                         time.taken <- end.time - start.time
                         print(time.taken)
                         print(paste("length rumors",length(getRumors(rud))))
                         rumsize=getRumors(rud)
                         qc=c()
                         qc$size=quantilediff(rumsize , obs$allrsize)
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
                         kc$size=kldiff(rumsize , obs$allrsize)
                         #kc$breadth=kldiff(rud$breadth , obs$allbreadth)
                         #kc$depth=kldiff(rud$depth , obs$alldepth+1)
                         kcrev=c()
                         #kcrev$size=kldiff(rumsize , obs$allrsize)
                         #kcrev$breadth=kldiff(obs$allbreadth , rud$breadth )
                         #kcrev$depth=kldiff(obs$alldepth+1,rud$depth)
                         if(full){
                             rud=list(rd=rud,t=time.taken)
                             save(rud,file=file.path(fold,"fullsim",paste0("res_sim_",i,".bin")))
                         }
                         print(list(kc=kc,qc=qc,kokc=kokc,kcrev=kcrev,kokcrev=kokcrev))
                         return(list(kc=kc,qc=qc,kokc=kokc,kcrev=kcrev,kokcrev=kokcrev))

                    },obs=observation,parameters=parameters,fold=fold,full=full)


    save(parameters,file=file.path(fold,"parameters.bin"))
    save(scores,file=file.path(fold,"scores.bin"))
}

mpi.close.Rslaves(dellog = T)
mpi.quit()

