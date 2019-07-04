### This simple script allows to try to run 50 simulation splitting those 50 simulation on 50 mpi slave
source("model.R")
library(Rmpi)

betamean=seq(0,1,length.out=10)  
u=seq(0,1,length.out=10)
mpi.spawn.Rslaves(nslaves=49)

mpi.bcast.cmd( id <- mpi.comm.rank() )
mpi.bcast.cmd( ns <- mpi.comm.size() )
mpi.bcast.cmd( host <- mpi.get.processor.name() )
mpi.bcast.cmd(source("model.R"))
mpi.bcast.cmd(source("model.R"))


normal_b=mpi.parLapply(1:50,function(i,u,betamean){
            nrt= lapply(betamean,function(b)cascades3D(N=10000,Nmax=400,IC=100,dtime=1,betadistrib=rnorm(10000,b),utility=u,R=10,repetition=500,lambda_c=1)[,c("size","U")])
},u=u,betamean=betamean)

save(normal_b,file="normal.bin")

mpi.close.Rslaves(dellog = FALSE)
mpi.quit()


