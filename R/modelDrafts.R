#' Simplemodel with only one message being retweeted
#' @param N : number of agents
#' @param K : max number of neighbours
#' @param betarange : range for parameter beta
#' @param utility : intrinsic utility of the message
#' @param repetition : number of time that the agents will check their timeline 
#' @param M : number of messageG
#'
#' @export simplemodel 


simplemodel <- function(agents=2:200,K=10,betadistrib=runif(0,1),jdistrib=runif(0,1),utility=runif(1),repetition=10,M=5,seed){

    N=length(agents)
    beta=betadistrib
    j=jdistrib
    tweets=rep(0,N)

    #listofneighbors = sapply(agents,function(a)sample(agents[-a],round(runif(1,1,K))))


    message=sample(agents,M)

    tweets[message]=1

    #print(paste("utility", utility,"rangebeta",paste0(range(beta),collapse="-")))

    for(i in repetition){
        nontweets=agents[!tweets[agents]]
        for(a in nontweets){
            randomneighbours=sample(agents[-a],round(runif(1,1,K)))
            k = length(randomneighbours) #total number of neighbours
            k_s=sum(tweets[randomneighbours]>= 1)      #number of neighbours that shared the message
            proba=exp(utility*beta[a]*k_s/k)/k
            if(is.nan(proba))proba=0
            #proba=exp(utility*beta[a]+k/k_s*j[a])/k
            if(proba > runif(1)){
                #print(paste("agent",a,"retweeted"))
                tweets[a]=max(tweets[randomneighbours]+1)
            }
        }
    }
    res=list()
    res$dist=sum(tweets)
    res$width=max(tweets)
    return(res)

}

#' Cascade number of retweet and length of cascad of differents tweets
#' @param N : number of agents
#' @param K : max number of neighbours (not used in this version)
#' @param M : number of messages
#' @param tl : maximum messages that agents can acess to
#' @param seed : initial number of messages 
#' @param seed : initial maximum number of tweeet that the `seeder` (first tweeter) can do 
#' @param betadistrib : vector os size N with distribution of beta for each agents
#' @param utility : vector of size M, the distribution of intrinsic utility for each messages
#' @param repetition : number of time that the agents will check their timeline 
#' @return a list with information on each messages  
#'
#' @export cascades 

cascades <- function(N=200,K=10,M=5,betadistrib=runif(200,0,1),jdistrib=runif(200,0,1),utility=seq(0,1,length.out=5),repetition=2,seed=100,tl=500){

    agents=1:N
    beta=betadistrib #should be a vector of float of size N
    j=jdistrib #should be a vector of float of size N

    tweets=matrix(0,nrow=N,ncol=M)
    messages=1:M
    names(messages)=1:M

    rownames(tweets)=agents
    colnames(tweets)=messages

    names(utility)=names(messages)

    #firstmover=sample(agents,seed) #choose the first agent sthat will tweet
   #seeder=sample(N,seed)   #choose message tweeted amoung all possible messages

    for(m in messages) #for all message
    {
        seeders=sample(N,seed)   #select number of firstmove
        tweets[seeders,m] = 1  #seeders tweet m
    }

    #print(paste("utility", utility,"rangebeta",paste0(range(beta),collapse="-")))

    for(r in 1:repetition){
        
        for(a in sample(agents)){

            #randomneighbours=sample(agents[-a],round(runif(1,1,K))) #to  generate random neighbors with K the max

            tweeted=messages[tweets[a,]>0] #messages that `a` has already tweeted

            subpossible=unique(cbind.data.frame(agent=sample((1:N)[-a],tl,replace=T),message=sample(M,tl,replace=T))) #if n and m are not big enough, this methods lead to sample twice the same tweet.

            subpossible = subpossible[apply(subpossible,1,function(param)tweets[param["agent"],param["message"]]>0),]#all messages indeed emitted 

            if(nrow(subpossible)>0){
                p=tapply(1:nrow(subpossible) , subpossible$message,length) #count number of messages of the same type

                probas =  exp(utility[names(p)] * beta[a]) #+ j[a] * p[i]  #selection function

                probas=probas/sum(probas) #normalize proba

                retweets=names(probas)[probas > runif(length(probas)) ] #new tweet with proba given by the score 
                #print(retweet)

                #increase cascade
                retweets=retweets[!(retweets %in% tweeted)] #eclude the tweet already tweeted
                for(t in retweets ){
                    ag = sample(subpossible$agent[ subpossible$message == t ],1) #select one of the people that tweeted `t`
                    tweets[a,t] = tweets[ag,t]+1 #get its cascage and increase it  
                } 

            }
            #print(probas)


        }
    }
    res=list()
    res$size=apply(tweets,2,function(t)sum(t>0))
    res$depth=apply(tweets,2,max)

    return(res)

}

