
args=commandArgs(trailingOnly = TRUE) #pass number of slave by comand line, should be #node - 1 as one node is used by the master

ns=args[1]
library(Rmpi)

betamean=seq(0,1,.05)  
u=seq(0,1,length.out=10)

mpi.spawn.Rslaves(nslaves=ns) 

mpi.bcast.cmd( id <- mpi.comm.rank() )
mpi.bcast.cmd( ns <- mpi.comm.size() )
mpi.bcast.cmd( host <- mpi.get.processor.name() )



normal_b=mpi.parSapply(betamean,function(b,u){
            print(paste("mean of beta=",b)) #output something for logging purpose
            nrt= sapply(1:100,function(i){Sys.sleep(.01);u/rnorm(b)}) ##somme dummy thing that should be run on every nodes
},u=u)

save(normal_b,file="normal.bin")

mpi.close.Rslaves(dellog = FALSE)
mpi.quit()


