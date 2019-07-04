source("../spread/R/tools.R")

###This script generate list of parameters and split thos list in two categorie to be able to slightly optimise resource usage. 
###The best wya to do it should use the function log(time)=1.66*log(Nmax*repetition)-16 to estimate the time used given the parameters

args=commandArgs(trailingOnly = TRUE) #pass number of slave by comand line, should be #node - 1 as one node is used by the master 

ns=args[1]#first argument is the number of slave
mainfold=args[2] #second argument = mainfolder to store the results
if(is.na(mainfold) | mainfold=="") mainfold="simulations"

print(paste0("Abc will be stored in mainfolder:",mainfold))
dir.create(mainfold)

### Define Prior Ranges ###
prior_N=c(15000,30000)
prior_R=c(1500,2000)
prior_beta=c(-100,-10,0,10,100)
prior_utility=c(-1,0,1)
prior_repetition=c(50,300)
prior_captl=c(1,600)
prior_mu_c=c(0,.035)
prior_IC=c(500,2000)
prior_Nmax=c(1000,10000)
prior_dtime=50:100
prior_stime=5:50

### Define Constants and Settings ###
nsim=768 #number of simulations
nsubfold=ns
full=F #Boolean to save or not the full simulations

lengthsmall=c()
lengthlarge=c()

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
    preparameters=list(
                    ### Define Prior Ranges ###,
                    N=sample(prior_N[1]:prior_N[2],nsim,replace=T),
                    R=sample(prior_R[1]:prior_R[2],nsim,replace=T),
                    repetition=sample(prior_repetition[1]:prior_repetition[2],nsim,replace=T),
                    captl=runif(nsim,prior_captl[1],prior_captl[2]),
                    betaDistrib=generalPartition(nsim,length(prior_beta)),#this contain the percentage for each class of agent 
                    utility=generalPartition(nsim,length(prior_utility)), #this contain the percentage for each class of utility 
                    mu_c=runif(nsim,prior_mu_c[1],prior_mu_c[2]),
                    IC=sample(prior_IC[1]:prior_IC[2],nsim,replace=T),
                    Nmax=runif(nsim,prior_Nmax[1],prior_Nmax[2]),
                    dtime=sample(prior_dtime,nsim,replace=T),
                    stime=sample(prior_stime,nsim,replace=T)
                    )

    #split and store it given estimate duration
    filteredparticules=preparameters$N > preparameters$Nmax & preparameters$IC < preparameters$R & preparameters$IC/preparameters$R > .2
    preparameters = lapply(preparameters,filter,filter=filteredparticules)
    preparameters$simtime=preparameters$Nmax * preparameters$repetition  
    parameters=lapply(preparameters,filter,filter= (preparameters$simtime > 1200000) )
    save(parameters,file=file.path(fold,"large_parameters.bin"))
    lengthlarge=c(lengthlarge,length(parameters$Nmax))
    parameters=lapply(preparameters,filter,filter= (preparameters$simtime <= 1200000) )
    save(parameters,file=file.path(fold,"small_parameters.bin"))
    lengthsmall=c(lengthsmall,length(parameters$Nmax))


}
print(mean(lengthsmall)*.35)
print(mean(lengthlarge)*.35)
print((mean(lengthlarge)+mean(lengthsmall)))


