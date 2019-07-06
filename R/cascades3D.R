#' Spread of cascades of retweet 
#' 
#' \code{cascades3D} Generate cascades of retweets from a initial population of rumors given 
#' 
#' main model of the package
#' @param N total number of agents
#' @param R initial number of rumors
#' @param captl maximum number of messages that agents can acess to 
#' @param IC initial number of cascades 
#' @param betadistrib vector of size N with the value of beta for each agents
#' @param utility vector of size R, with the value of intrinsic utility for each rumors
#' @param time number of time that the agents will check their timeline  (length of the simulation)
#' @param stime the oldest tweet (in timestep) that an agent can access from an active cascade (if stime = 1 the agents can only access retweet or tweet done during the previous time step. If -1 there is no limit.
#' @param dtime maximum number of timestep during which a cascades can survives. if =< 0 cascades are kept during all the simulation.
#' @param mu_r rate of apparition of new rumors
#' @param mu_c rate of apparition of new cascades
#' @param summary if TRUE, return a data.frame with summary statistics about the cascade, else return a list of lists with a list of all cascades, a list of all rumors, a list of the results, and much more
#' @param Nmax if NULL or Nmax == 1, N agents will have the opportunity to retweet at each time step, ifNmax > 1, Nmax randomly sampled agent will have the possibility to retweet, if Nmax < 1, Nmax represnet the proportion of agents that will be able to retweet
#' @param memory if TRUE agents cannot participate to the same cascades twice (ie retweet twice the same tweet). Real way twitter is working.
#' @param log print the time step and some information if TRUE, nothing if FALSE
#' @return depends on the value of the BOOLEAN summary 
#' @examples
#' simplerun=cascades3D()
#' par(mfrow=c(1,2))
#' plotCCFD(simplerun$size,main="Cascade Size")
#' plotCCFD(simplerun$depth,main="Cascade Depth")
#' @export cascades3D 
cascades3D <- function(N=200,R=10,betadistrib=rep(1,200),utility=seq(0,1,length.out=10),time=2,captl=-1,IC=50,summary=T,Nmax=NULL,mu_c=0,mu_r=0,dtime=0,stime=-1,memory=T,log=T,metrics=c("size","depth","breadth"),debug=F){
        
    

    stopifnot(length(betadistrib)== N,length(utility)== R)
    if(summary && !memory && (sum(metrics != "size")>0) )stop("if memory is false depth and breadth will generate infinite loop when trying to get summary statistics")

    if(is.null(Nmax))Nmax=1
    if(Nmax<=1)Nmax=Nmax*N


    if(memory) #if we keep memory of previous tweets then we create a list that stores the cascades in which every agents have already participated
        agents=lapply(1:N,function(i)c()) 
    else  #if not, we keep a simple vector with number of agents
        agents=1:N
    beta=betadistrib
    names(agents)=beta #the list (or the vector) is named given the $\beta$-value of each agent
    rumors=1:R

    round=0

    newc=generateCascades(IC,N,R,round) #generate a new pool of cascade 
    cascades=newc$cascades
    cascades_id=newc$cascades_id
    if(memory){
        seeders=unique(cascades_id$seeder)
        for(s in seeders){
            id=which(cascades_id$seeder == s)
            agents[[s]]=c(agents[[s]],id)
        }
    }  

    activepool=1:nrow(cascades_id)

    for(r in 1:time){


        ##get all possible tweets: if done here all agent can see what the other agent have done previously during the same time step
        #available_t =lapply(cascades[activepool],function(c)c[,2]) #get the node of all cascade that can be retweeted
        if(stime<0)
            laststep=0
        else
            laststep=max(0,r-1-stime)
        if(log)print(paste("timestep",r,"round",round,"last round",laststep))
        available_t =lapply(cascades[activepool],function(ca)ca[ca[,3]>=laststep,2]) #get the node of all cascade that can be retweeted
        possible=unlist(available_t,use.names=F)        #create a 1D vector of all possible
        actives_id=seq_along(cascades)[activepool]      #among the list of cascades ID, get the one that are still actives
        #print(paste("twetosome size:",paste0(lengths(available_t),collapse=","),"cascadosome",nrow(cascades_id),"activisome",paste0(actives_id,collapse=","))  )
        possible_cid=rep(actives_id,lengths(available_t)) #it may be possible to directly sample using the lengths of the cascade as a proba in sample(x,n,proba), and then randomly select among the leaf of the cascade selected. But in that cas it's hard to know how many tweets agents can select using the parameter tl
        possible_t=length(possible)
        tl=captl
        #print(paste0("possibility:",possible_t," reall:",length(possible)," tl:",tl))
        if(captl<0)tl=possible_t
        else if(captl<1)tl=tl*possible_t
        else if(captl > possible_t) tl=possible_t
        #print(paste0("possibility:",possible_t," reall:",length(possible)," tl:",tl))

        if(debug)newretwett=c()
        ####
        for(a in sample.int(N,Nmax)){
            if(runif(1)>mu_c){
            subsample=sample.int(possible_t,tl) #select a subsample of the size tl if possible among all possible retweet
            tsource=possible[subsample]                #get the agent id of the source of the possible tweet
            tsource_cid= possible_cid[subsample]       #get the id of the cascade associated with this tweet
            tsource_rid=cascades_id$rumor[tsource_cid] #get the id of the rumor associated with this cascade

            u_subsample=utility[tsource_rid]           #get the utilities of the rumors involved in those cascade
            ##get all possible tweets: if done here all agent can see what the other agent have done previously during the same time step
            #if done before the for loop it means that all agents won have the information of wha thave been done by the other agent before the end of the current time step
            #available_t =sapply(cascades,function(c)c[,2]) #get the node of all cascade that can be retweeted
            #possible=unlist(available_t,use.names=F)        #create a 1D vector of all possible
            #possible_cid=rep(1:nrow(cascades_id),lengths(available_t)) #create a vector with the indice of the corresponding cascades
            ####


            ##compute proba
            p_i=probaSelection(u_subsample,beta[a])

            #print(paste0("possibility:",possible_t," selected:",length(p_i)))
            rnd=runif(length(p_i))
            #print(rbind(rnd,p_i))
            selected_rt=(p_i >rnd ) # & !(tsource_cid %in% tweeted)  #probability & not already participate 
            round=round+1
            npossibletweet=sum(selected_rt)
            if(debug)print(paste("nrt:",npossibletweet))
            if(debug)newretwett=c(newretwett,npossibletweet)


            if(npossibletweet>0){# the agent may retweet something
                ### The following ifs statement are here to avoid some useless and "potentially"  costly computations
                if(memory){
                    tweeted=agents[[a]]     #get id of cascade in which agent already participated
                    #print(paste("agents",a,"gonna selected:",paste0(tsource_cid[selected_rt],collapse=",")))
                    #print(paste("agents",a,"already tweeted:",paste0(tweeted,collapse=",")))
                    if(length(tweeted)>0){ #if we already tweet something then it's worth looking at what
                        not_already_tweeted = !(tsource_cid %in% tweeted)   #check which among possible tweet haven't been tweeted yet
                        #print(tsource_cid)
                        #print(not_already_tweeted)
                        if(sum(not_already_tweeted)>0){# the agent haven't already retweeted everything among the possible tweets he sees
                            newselection=selected_rt & not_already_tweeted
                            #print(newselection)
                            if(sum(newselection)<npossibletweet){
                                selected_rt=newselection
                                #print(paste("finally agents",a,"gonna select:",paste0(tsource_cid[selected_rt],collapse=",")))
                                npossibletweet=sum(selected_rt)
                            }
                        }
                        else{ #agent have already participated to all the available cascades.
                            npossibletweet=0
                        }

                    }
                }
                if(npossibletweet>0){#we check again if the agent will indeed retweet something
                    new_retweet_cid=tsource_cid[selected_rt]
                    new_retweet=tsource[selected_rt]

                    if(sum(duplicated(new_retweet_cid))>0){ #if there are some duplicate a few checks have to be done. They may be costly that's why we avoid them  as much as possible
                        nodouble=removeDuplicate(new_retweet_cid) #we can only have 1 retweet per cascade thus in the case that we picked up two tweet from the same cascade, we randomly choose one of those tweet to retweet. This can be costly in theory as it goes through all the list new_retweet_cid, but in practice this list should stay relatively small 
                        new_retweet_cid=new_retweet_cid[nodouble]
                        new_retweet=new_retweet[nodouble]
                    }

                    #print(paste("timestep",r,"agent",a,"will retweet",paste0(new_retweet,collapse=","),"from cascades",paste0(new_retweet_cid,collapse=","),"at round",round))

                    if(memory)agents[[a]]=c(agents[[a]],unique(new_retweet_cid)) #add this cascade among the cascade retweeted by the agent
                    cascades[new_retweet_cid]=lapply(1:length(new_retweet_cid),function(c)rbind(cascades[[new_retweet_cid[c]]],c(new_retweet[c],a,round)))
                }

            }
            }
        }
        if(debug)print(paste("mean number of retweet per user: ",mean(newretwett)))
        prevnc=nrow(cascades_id) #we count the number of existing cascade, will be used to give an ID to the new one
        nc=sum(runif(Nmax)<mu_c)
        #nr=rpois(1,mu_r)
        if(nc>0){
            newc=generateCascades(nc,N,R,round)
            if(memory){
                seeders=unique(newc$cascades_id$seeder)
                for(s in seeders){ #here same agent can introduce again the exact same rumor.
                    id=which(newc$cascades_id$seeder == s)
                    agents[[s]]=c(agents[[s]],id+prevnc)
                }
            }
            cascades=c(cascades,newc$cascades)
            cascades_id=rbind(cascades_id,newc$cascades_id)
        }

        if(log)print(paste("round", round, "ncascades:",length(cascades),"actives:",length(activepool[activepool])))

        ##Next lines handle if cascade can be "deactivate" when they haven't been retweet during a certain number of timestep dtime
        #the `sapply` here may be better done by not going through the _already desactivated_ cascades. Maybe something like: `seq_along(cascades)[activisome]`, as done before
        if(dtime>0)
            activepool=sapply(cascades,function(ca)max(ca[,3])>(Nmax*(r-dtime))) #keep only the cascade active during the last round of retweet
        else
            activepool=1:(length(cascades))
        stopifnot(length(cascades) == nrow(cascades_id))

        #Take care of extinction case
        if(length(activepool[activepool])==0){
            if(log)print("nothing left to retweet");
            if(summary)  cbind.data.frame(U=utility[cascades_id$rumor],cascades_id,getSummaryCascade(cascades,metrics=metrics)) #count depth and size and add that to general dataframe
            else list(agents=agents,summary=cbind.data.frame(U=utility[cascades_id$rumor],cascades_id),cascades=cascades)
        }

    }

    if(summary)  cbind.data.frame(U=utility[cascades_id$rumor],cascades_id,getSummaryCascade(cascades,metrics=metrics)) #count depth and size and add that to general dataframe
    else list(agents=agents,summary=cbind.data.frame(U=utility[cascades_id$rumor],cascades_id),cascades=cascades)
}

