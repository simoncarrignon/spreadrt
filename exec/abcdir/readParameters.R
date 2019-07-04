source("../spread/R/cascades3D.R")
source("../spread/R/tools.R")
source("../spread/R/ABC_analyse.R")

args=commandArgs(trailingOnly = TRUE) #pass number of slave by comand line, should be #node - 1 as one node is used by the master 

paramfold=args[1] #first argument = paramfold where the read parameter are stored 
mainfold=args[2] #second argument = mainfolder to store the results
timecat=args[3] #length of the simulation (categorie) 
ns=args[4] #number of slave used 
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

### Define Constants and Settings ###
full=F #Boolean to save or not the full simulations


load("observations.bin")

fi=0
load(file.path(paramfold,paste0(timecat,"_parameters.bin")))
parameters$simtime=NULL
nsim=length(parameters$Nmax)

fold=file.path(mainfold,paste0(mainfold,fi))
while(file.exists(fold)){
    fi=fi+1
    fold=file.path(mainfold,paste0(mainfold,fi))
}
dir.create(fold)
if(full)dir.create(file.path(fold,"fullsim"))

print(paste(timecat,paramfold))

scores=mpi.applyLB(1:nsim,function(i,parameters,obs,fold,full){
                   start.time <- Sys.time()
                   param=lapply(parameters,filter,filter=i)
                   print(i)
                   print(param)
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

mpi.close.Rslaves(dellog = T)
mpi.quit()

