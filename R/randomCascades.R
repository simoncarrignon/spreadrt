#' randomCascades: A slightly (almost not) modified version of the original model by Alex & Damien
#' Original model is in tools/FNM.R
#' @return dataframe that allow to compare easily to the result of the model of the cascade
#' @export randomCascades 
randomCascades <- function(Nmin=100000,Nmax=50000,t_steps=50,y_max=1,runs=1,mu=0.00001,tau = 1,alpha=0,conformity=F,beta=-2,topfive=F,TF=5,C=.1,alberto=F){

    t<-seq(1,t_steps)
    N_vec=round(Nmin*exp(alpha*t))

    T <- as.integer(tau + t_steps) # total time_steps

    N=as.integer(N_vec[1])#inital N value

    old_pop <- seq(1,N) #inialise population
    base = N+1

    corpus <- list()
    vocab<- list()
    Z <-matrix(nrow=y_max,ncol=t_steps-1)

    #save all the values in one popualtion
    total_pop <- vector(mode='integer',length=N*t_steps)

    #beginning of t loop 
    for (t in 1:T){

        rnd = runif(N, min = 0, max = 1)

        ##### the assingment was set to pop, i've changed it to old_pop
        select_vec<-old_pop 

        pop<-old_pop # start new empty population


        if(alberto){
            pop=originalTopfive(pop,TF,C,rnd,mu,base)
        }
        else{ 
            if(conformity){ #introduction of a simple conformity bias 
                #freq=table(pop)
                #rnd=rnd/(freq[as.character(pop)])
                Bias <- conformityBias(old_pop,beta)
                pop[rnd>=mu] <- sample(old_pop, length(pop[rnd>=mu]), replace = TRUE,prob = Bias)
            }
            else{
                pop[rnd>=mu] <- sample(old_pop, length(pop[rnd>=mu]), replace = TRUE)
            }
            if(topfive){
                candidates=topfive(old_pop,TF,C)[runif(length(old_pop)*C)>=mu]
                pop[sample.int(length(old_pop),length(candidates))]= candidates
            }
        }


        pop[rnd<mu] <- seq(base,base+length(pop[rnd<mu])-1)


        base <- base + length(pop[rnd<mu]) # next trait to enter the population
        if (t > tau){
            N <- as.integer(N_vec[t-tau])
            s<-t-tau-1
            first<-(1+(s*N))
            last<-((s*N)+N)
            total_pop[first:last]<-pop
        }

        #rank_old <- rank_now
        old_pop <- pop # current pop becomes old_pop
        #end of t loop
    }

    MC <- sort(table(total_pop),decreasing = TRUE)
    df<-data.frame(rank=rank(MC,ties.method = "first"),MC)

    colnames(df)=c("rank","U","size")
    df$U=as.numeric(df$U)
    df
}

#' @export conformityBias 
conformityBias <- function(pop,beta) {
  
    #frequency Pk of each trait
    Pk<-(table(pop))
    Pk<-Pk/sum(Pk)
    Pk<-data.frame(Pk)

    #the probability of each entry in the pop
    probs <- pop
    probs[] <- Pk$Freq[match(unlist(pop), Pk$pop)]
    top=T

    if(top){
        #apply the exponential 
        exp<-exp(beta*probs)
        probs[probs>0.001]<-exp/sum(exp)
    }
    else{
        #apply the exponential to everyone 
        exp<-exp(beta*probs)
        probs<-exp/sum(exp)
    }
  
  
  return(probs)
}

#' @export originalTopfive 
originalTopfive <- function(pop,TF,C,rnd,mu,base){
        #the most popular traits in the population
    old_pop=pop
        topTraits <- sort(table(pop),decreasing=TRUE)[1:TF]
        topTraits <- as.data.frame(topTraits)$pop


        pop[rnd>=mu & rnd<(1-C)*(1-mu)] <- sample(old_pop, length(pop[rnd>=mu]), replace = TRUE)    
        pop[rnd<mu] <- seq(base,base+length(pop[rnd<mu])-1)
        
        
        #the section of the population that has information about the most popular trair
        popConform <- pop[rnd > (1-C)*(1-mu)]
        
        # only copy if you do not have a top trait
        Ncopy<-length(popConform[popConform != topTraits[!is.na(topTraits)]]) #SC: i had to add this test because some time topTraits include a NA. It looks strange and it should not be possible.
        popConform[popConform != topTraits[!is.na(topTraits)]] <- sample(old_pop, Ncopy, replace = TRUE) 

        
        pop[rnd > (1-C)*(1-mu)]<-popConform
        pop
}



#' @export topfive 
topfive <- function(pop,TF,C){ 

    tpop= table(pop)
    topfive=sort(tpop,decreasing=T)[1:TF]
    topfive=as.data.frame(topfive)$pop
    cpop=sample(topfive,length(pop)*C,replace=T)
    return(cpop)

}


