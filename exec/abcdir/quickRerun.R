source("../spread/R/cascades3D.R")
source("../spread/R/randomCascades.R")
source("../spread/R/tools.R")
source("../spread/R/ABC_analyse.R")

#args=commandArgs(trailingOnly = TRUE) #pass number of slave by comand line, should be #node - 1 as one node is used by the master 
#
#ns=args[1]#first argument is the number of slave
#mainfold=args[2] #second argument = mainfolder to store the results
#if(is.na(mainfold) | mainfold=="") mainfold="simulations"
#
print("rerun nonneutrel model be stored in mainfolder:")


library(Rmpi)
mpi.spawn.Rslaves(nslaves=2)

mpi.bcast.cmd( id <- mpi.comm.rank() )
mpi.bcast.cmd( ns <- mpi.comm.size() )
mpi.bcast.cmd( host <- mpi.get.processor.name() )
mpi.bcast.cmd(source("../spread/R/randomCascades.R"))
mpi.bcast.cmd(source("../spread/R/cascades3D.R"))
mpi.bcast.cmd(source("../spread/R/tools.R"))
mpi.bcast.cmd(source("../spread/R/ABC_analyse.R"))
mpi.bcast.cmd(mainfold)

    size=list()
    size$true=trueru$size
    size$false=falseru$size
    bins=10^(seq(0,5.5,.1))
    tfcols=c("green","red")
    names(tfcols)=c("true","false")
    tf=c("true","false")
    names(tf)=c("true","false")


    load("nonneutral.posteriors.bin")
    scores=mpi.applyLB(1:nsim,function(i,parameters,obs,fold,full){
                       reruncontent=lapply(tf,function(t)lapply(c(3,5),function(i)getRumors(cascades3D.list(unlist(nonneutral.posteriors$qc[[t]][i,1:17]),metrics="size"))))
                       save(file="reruncontent",reruncontent)

mpi.close.Rslaves(dellog = T)
mpi.quit()
