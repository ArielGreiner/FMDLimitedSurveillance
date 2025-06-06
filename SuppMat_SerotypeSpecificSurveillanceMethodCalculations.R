####CODE TO ASSESS THE SURVEILLANCE METHODS FOR SEROTYPE A AND O AND GENERATE THE FIGURES IN THE SUPPLEMENTARY MATERIAL

#probably don't need all of these...
library(scales)
library(deSolve)
library(geometry)
library(fields)
#library(optparse)
library(maps)
library(rgdal)
library(ggplot2)
library(viridis)
library(dplyr)
library(plyr)
library(igraph)
library(purrr)
library(tidyverse)

####SET WORKING DIRECTORY####
setwd('/Users/akg6325/Dropbox/Github') #replace with your own directory!

####LOAD IN USEFUL FUNCTIONS####

#functions written by José 

### Function that, given a network object (igraph), calculates node level statistics
### Works for networks of any size
nodeStatsDW<-function(x){
  eig<-eigen_centrality(x,directed = T,weights = E(x)$weight)$vector
  inDegree<-degree(x,mode = "in")
  outDegree<-degree(x,mode = "out")
  Degree<-degree(x,mode="total") #added this one because i want to use it :)
  inStrength<-strength(x,mode="in",weights = E(x)$weight)
  outStrength<-strength(x,mode="out",weights = E(x)$weight)
  betweenness<-betweenness(x, directed=T, weights=E(x)$weight,nobigint = TRUE, normalized = FALSE) #log(1/c_ij to approx distance?)
  inCoreness<-coreness(x,mode="in")
  outCoreness<-coreness(x,mode="out")
  centrality <- data.frame(Node = as.numeric(rownames(as.data.frame(inDegree))),
                           OutDegree    = outDegree,
                           InDegree    = inDegree,
                           Degree = Degree,
                           OutStrength = outStrength,
                           InStrength   = inStrength,
                           OutCoreness = outCoreness,
                           InCoreness = inCoreness,
                           Betweenness  = betweenness,
                           EigenCent = eig
  )
  centrality
}

### Function to create directed weighted networks. Need a file with 3 columns"
### Source, Target and Weight
### This works for networks of any size
networkCreationDW<-function(x){
  el1 <-as.data.frame(x) #Verify that it's a data.frame
  #  el1<-subset(el1,source_epiunit_id!=destination_epiunit_id) #Avoid loops, self-connections
  names(el1)<-c("i","j","w")
  el1<-subset(el1,i!=j) #Avoid loops, self-connections
  #  el1<-subset(el1,Source!=Target) #Avoid loops, self-connections
  ell3<-as.matrix(el1)
  ell3[,1]=as.character(ell3[,1]) #Because the vertex IDs in this dataset are numbers, we make sure igraph knows these should be treated as characters. Otherwise, it'll create problems (see page on data import)
  ell3[,2]=as.character(ell3[,2])
  g=graph.edgelist(ell3[,1:2],directed = TRUE) #Create the directed network
  E(g)$weight=as.numeric(ell3[,3]) #weights using the number of cattle made
  g<-simplify(g) #Verify that there are not multiple connection or remaining loops
}

#Function that works with a network in igraph format
#creates a file with some measures
netsRealStats<-function(g){
  compsG<-clusters(g)
  x<-induced.subgraph(g,compsG$membership==which.max(compsG$csize))
  pi<-transitivity(x,type="weighted",weights = E(x)$weight)
  is.na(pi)<-sapply(pi,is.infinite)
  return(c(vcount(g),ecount(g),vcount(x),ecount(x),mean(pi,na.rm=T),
           eigen_centrality(x,directed = T,weights = E(x)$weight)$value,
           reciprocity(x),
           mean_distance(x, directed = TRUE),
           assortativity(x, types1 = graph.strength(x,mode = "in"),
                         types2 = graph.strength(x,mode="out"),directed = T),
           diameter(x,directed = TRUE,weights =E(x)$weight)))
  #  statsGlobalReal<-reci
  #  write.table(reci,file="~/Projects/FMD-Surveillance/raw-data/FrequencyRelated/RealCompleteNetworkFrequency.csv",append = TRUE,
  #              sep = ",",col.names=FALSE,row.names = F)
  #  statsGlobalReal
}


####LOAD IN OUTBREAK AND NETWORK DATA####
#outbreak data
outbreakepiunits <- readRDS(file = "FMDLimitedSurveillance/Data/outbreakepiunits.rds") 
outbreakData <- readRDS(file = "FMDLimitedSurveillance/Data/outbreakData_final.rds") 
#note that the position_x and position_y columns were randomized for data sensitivity reasons
#chose random values for each vill_ID that ranged from 0 -> the range of x or y values respectively
#newlocations <- data.frame(vill_ID = seq(1,55193,1), x = runif(55193, min = 0, max = range_x), y = runif(55193, min = 0, max = range_y))

t_outbreaks <- readRDS(file = paste0("FMDLimitedSurveillance/Data/t_outbreaks.rds"))
tplus1stats_full_fxd <- readRDS(file = paste0("FMDLimitedSurveillance/Data/tplus1stats_full_fxd.rds"))

#2 month network data
#created from an edge list using the networkCreationDW code above
networks_2monthsep <- readRDS(file = "FMDLimitedSurveillance/Data/networks_2monthsep.rds") 
#created from an edge list using the nodeStatsDW code above
networkmetrics_2monthsep <- readRDS(file = "FMDLimitedSurveillance/Data/networkmetrics_2monthsep.rds") 

#location data
#note that the position_N and position_E columns were randomized for data sensitivity reasons
#chose random values for each vill_ID that ranged from 0 -> the range of E (corresponds with x) or N (corresponds with y) values respectively
#newlocations <- data.frame(vill_ID = seq(1,55193,1), x = runif(55193, min = 0, max = range_x), y = runif(55193, min = 0, max = range_y))
locationsmasterlist <- readRDS(file = "FMDLimitedSurveillance/Data/locationsmasterlistabr_final.rds") 

#all time network data
epiunit_networkmetrics <- readRDS(file = "FMDLimitedSurveillance/Data/epiunit_networkmetrics.rds")



####SPATIAL PROXIMITY METHOD####
#need to load in: outbreakData, outbreakepiunits_A, outbreakepiunits_O, locationsmasterlist

#run a for loop where, for every month, the radius slowly increases until it exceeds x%*(totalepiunits) and then it records data from that and continues until it hits >35% of all epiunits

steps <- 500
pct <- c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.5,0.75,1)

numserotypes <- 2

for(n in 1:numserotypes){
  if(n == 1){
    outbreakepiunits <- outbreakepiunits_A
    serotype <- "A"
  }
  
  if(n == 2){
    outbreakepiunits <- outbreakepiunits_O
    serotype <- "O"
  }
  
  
  euclideandist_varr_outneighbs <- data.frame(twomthntwk = rep(seq(1,num_2_ntwks,1), each = length(pct)), pctaway = rep(pct, num_2_ntwks), pctntwk = NA, avgsize_peroutbk = NA, pcttplus1outbks = NA, radii = NA)
  
  proximity_varr_tplus1epiunits <- list()
  
  
  for(i in 1:num_2_ntwks){
    print(paste("i = ", i))
    t_outbreakdata <- locationsmasterlist[(locationsmasterlist$vill_ID %in% outbreakepiunits[[i]]),]
    if(dim(t_outbreakdata)[1] == 0){next}
    counter <- 0
    for(j in 1:steps){
      stepsize <- j/50 #before steps = 20, and stepsize was j/4, so hopefully this will be way too small in a good way
      for(k in 1:length(outbreakepiunits[[i]])){
        outneighbourhoods <- list()
        outneighbourhoods[[k]] <- locationsmasterlist$vill_ID[(locationsmasterlist$position_N > (t_outbreakdata$position_N[k]-stepsize)) & (locationsmasterlist$position_N < (t_outbreakdata$position_N[k]+stepsize)) & (locationsmasterlist$position_E > (t_outbreakdata$position_E[k]-stepsize)) & (locationsmasterlist$position_E < (t_outbreakdata$position_E[k]+stepsize))]
        if(k == 1){epiunits_outneighbourhoods <- outneighbourhoods[[1]]}
        if(k > 1){
          epiunits_outneighbourhoods <- append(epiunits_outneighbourhoods, outneighbourhoods[[k]][which(!(outneighbourhoods[[k]] %in% epiunits_outneighbourhoods))])
        }
      }
      for(l in 1:length(pct)){
        if(l > counter){ #go through only the radii that haven't been filled in yet
          if((length(epiunits_outneighbourhoods)) >= (numtotalepiunits*pct[l])){ #if see greater than or = the number of epiunits at that percentage, fill in the information
            
            euclideandist_varr_outneighbs$avgsize_peroutbk[euclideandist_varr_outneighbs$twomthntwk == i & euclideandist_varr_outneighbs$pctaway == pct[l]] <- sum(lengths(outneighbourhoods))/length(outbreakepiunits[[i]])
            
            euclideandist_varr_outneighbs$pctntwk[euclideandist_varr_outneighbs$twomthntwk == i & euclideandist_varr_outneighbs$pctaway == pct[l]] <- length(epiunits_outneighbourhoods)/numtotalepiunits
            
            euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$twomthntwk == i & euclideandist_varr_outneighbs$pctaway == pct[l]]  <- length(epiunits_outneighbourhoods[(epiunits_outneighbourhoods %in% outbreakepiunits[[(i+1)]])])/length(outbreakepiunits[[(i+1)]])
            
            euclideandist_varr_outneighbs$radii[euclideandist_varr_outneighbs$twomthntwk == i & euclideandist_varr_outneighbs$pctaway == pct[l]] <- stepsize
            
            proximity_varr_tplus1epiunits[[((i-1)*(length(pct)) + l)]] <- epiunits_outneighbourhoods[(epiunits_outneighbourhoods %in% outbreakepiunits[[(i+1)]])]
            
            counter <- l #if go into this if statement, mean you're > goal % age for pct[l] so don't want to go in there again
            print(paste("counter =", counter, "stepsize = ", stepsize))
            if(counter == length(pct)){break}
          }
        }
      }
      if(counter == length(pct)){break} #don't need to keep going through the steps, this only breaks you out of the forloop it's in..not the whole nested set
    }
    counter <- NULL
    outneighbourhoods <- epiunits_outneighbourhoods <- NULL
  }
  
  #saveRDS(euclideandist_varr_outneighbs, file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/SerotypeSpecific/Serotype",serotype,"/2months_euclidean_varr_distoutneighbourhoodmetrics_11.26.2023.rds")) 
  #saveRDS(proximity_varr_tplus1epiunits, file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/SerotypeSpecific/Serotype",serotype,"/2months_proximityvarr_tplus1epiunits_11.26.2023.rds")) 
  
  serotype <- outbreakepiunits <- NULL
}

####NETWORK CONNECTIVITY METHOD####
#need to load in networkmetrics_2monthsep, outbreakepiunits_A, outbreakepiunits_O, epiunit_networkmetrics, outbreakData

#divide by the number of outbreaks observed
p_perc = c(0.05,0.1,0.15,0.20,0.25,0.3,0.35)
metrix <- c("Degree")
numtotalepiunits <- 54096

numserotypes <- 2
for(n in 1:numserotypes){
  
  if(n == 1){
    outbreakepiunits <- outbreakepiunits_A
    serotype <- "A"
  }
  
  if(n == 2){
    outbreakepiunits <- outbreakepiunits_O
    serotype <- "O"
  }
  
  pperc_outbreaks_full <- data.frame(metric = rep(rep(metrix, each = length(p_perc)),(num_2_ntwks+1)),P_perc = rep(rep(p_perc,length(metrix)), (num_2_ntwks+1)), twomonthpd = rep(c(seq(1,num_2_ntwks,1), "alltime"), each = (length(p_perc)*length(metrix))), Actual_peroutbreaks = NA)
  
  
  for(j in 1:num_2_ntwks){
    topDegreeepiunits <- list()
    networkmetrics_2monthsep_indiv <- networkmetrics_2monthsep[[j]] 
    onemonthoutbreaks <- outbreakepiunits[[(j+1)]]
    
    for(i in 1:length(p_perc)){
      #which are the top P% of nodes for each metric?
      topDegreeepiunits[[i]] <- networkmetrics_2monthsep_indiv$Node[order(networkmetrics_2monthsep_indiv$Degree, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)] 
      
      pperc_outbreaks_full$Actual_peroutbreaks[pperc_outbreaks_full$metric == "Degree" & pperc_outbreaks_full$P_perc == p_perc[i] & pperc_outbreaks_full$twomonthpd == j] <- length(which(topDegreeepiunits[[i]] %in% onemonthoutbreaks))/length(onemonthoutbreaks) 
    }
    #erase all of these, just to make sure everything is working
    topDegreeepiunits <- NULL
  }
  
  
  #want to compare that with the full network/full timeseries

  #which are the top P% of nodes for each metric?
  #p_perc = c(0.1,0.25,0.5,0.75)
  topDegreeepiunits <- list()
  for(i in 1:length(p_perc)){
    topDegreeepiunits[[i]] <- epiunit_networkmetrics$Node[order(epiunit_networkmetrics$Degree, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)]
  }
  
  
  #what percentage of the nodes in the top P% for each metric, experience outbreaks?
  
  #make a vector containing the epiunit # of all of the epiunits that got infected
  #outbreakData probably/might have duplicates in vill_ID
  outbreakepiunits_nontemp <- unique(unlist(outbreakepiunits)) #unique(outbreakData$vill_ID) <- non serotype specific version
  
  #make that table, finally
  #Actual = any epiunits that had an outbreak at ANY point, since infections don't disappear in this system
  #divide by the number of outbreaks observed
  metrix <- c("Degree")
  pperc_outbreaks <- data.frame(metric = rep(metrix, each = length(p_perc)),P_perc = rep(p_perc,length(metrix)), Actual_peroutbreaks = NA)
  
  for(i in 1:length(p_perc)){
    pperc_outbreaks$Actual_peroutbreaks[pperc_outbreaks$metric == "Degree" & pperc_outbreaks$P_perc == p_perc[i]] <- length(which(topDegreeepiunits[[i]] %in% outbreakepiunits_nontemp))/length(outbreakepiunits_nontemp) #OLD: length(which(outbreakepiunits_nontemp %in% topDegreeepiunits[[i]]))/length(outbreakepiunits_nontemp)
    
  }
  
  #add it to the temporal networks dataframes
  pperc_outbreaks_full$Actual_peroutbreaks[pperc_outbreaks_full$twomonthpd == "alltime"] <- pperc_outbreaks$Actual_peroutbreaks
  
  #saveRDS(pperc_outbreaks_full, file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/SerotypeSpecific/Serotype",serotype,"/2monthnetworks_degreeonlyriskfactor_fullnetworkvalues_11.26.2023.rds")) 
  
  serotype <- outbreakepiunits <- pperc_outbreaks_avg <- pperc_outbreaks_full <- NULL
}

####NETWORK PROXIMITY METHOD####
#need networks_2monthsep, t_outbreaks_a/_o (because only looking at the t outbreaks that act as sources in the network), tplus1_stats_a/_o (because only going to find outbreaks that have network destination data), outbreakepiunits

numtotalepiunits <- 54096 #from: dim(locations)[1] see 9.21.2023 chunk above
pct_smol <- c(0.05,0.1,0.15,0.2,0.25,0.3,0.35)
maxorder = 20

numserotypes <- 2
for(n in 1:numserotypes){
  
  if(n == 1){
    outbreakepiunits <- outbreakepiunits_A
    tplus1stats_full_fxd <- tplus1_stats_a
    t_outbreaks <- t_outbreaks_a
    serotype <- "A"
  }
  
  if(n == 2){
    outbreakepiunits <- outbreakepiunits_O
    tplus1stats_full_fxd <- tplus1_stats_o
    t_outbreaks <- t_outbreaks_o
    serotype <- "O"
  }
  
  wctracing_outneighbs <- data.frame(twomthntwk = rep(seq(1,num_2_ntwks,1), each = length(pct_smol)), pctaway = rep(pct_smol, num_2_ntwks), pctntwk = NA, pcttplus1outbks = NA, edgerank = NA, order = NA)
  wcontacttracing_tplus1epiunits <- list()
  
  
  start_time <- Sys.time()
  for(i in 1:num_2_ntwks){
    print(paste("i = ", i))
    tplus1outbreaks <- tplus1stats_full_fxd[[i]]$epiunit[tplus1stats_full_fxd[[i]]$inf_status == 1 & tplus1stats_full_fxd[[i]]$tplus1 == (i+1)]
    counter <- 0
    totalnodes <- NULL
    srchdnodes <- NULL
    for(j in 1:maxorder){
      if(j==1){#searching the edges of the t outbreak nodes
        sourcenodes <- t_outbreaks[[i]]
        srchdnodes <- sourcenodes #want to keep track of the nodes whose outgoing links have been searched
      }
      if(j > 1){ #searching the epiunits in the outneighbourhoods of the previous steps sourcenodes, but only looking at nodes that haven't been searched before 
        sourcenodes <- newsources[which(!(newsources %in% srchdnodes))]
        #namepi <- names(epiunits_outneighbourhoods)
        #sourcenodes <- as.integer(namepi[which(!(namepi %in% as.character(srchdnodes)))])
        srchdnodes <- append(srchdnodes, sourcenodes)
        newsources <- NA
        #epiunits_outneighbourhoods <- NA
        #outneighbourhoods <- NA
      }
      if(length(sourcenodes) == 0){
        sourcenodes <- mostedges <- rkd_epiunits <- NA
        break
      } #if no more source nodes, go to the next loop
      #rank the epiunits each of the source nodes are linked to, based on the weight of their edge 
      #first need to make a dataframe to search through later
      mostedges <- max(degree(networks_2monthsep[[i]],v=as.character(sourcenodes),mode = "out")) #give each source node the max # of possible candidate node spots
      if(mostedges == 0){ #if the largest degree of any of the sourcenodes is 0...there will be nothing to search
        sourcenodes <- mostedges <- rkd_epiunits <- NA
        break
      } 
      rkd_epiunits <- data.frame(source = rep(sourcenodes, each = mostedges), rank = rep(seq(1,mostedges,1), length(sourcenodes)), candidates = NA)
      for(a in 1:length(sourcenodes)){
        edgz <- E(networks_2monthsep[[i]])[.from(as.character(sourcenodes[a]))] #doesn't include self loops for some reason
        epi_rank <- order(edgz$weight, decreasing = TRUE) #for ties, it goes in the order that they were listed (not sure how that's decided)
        if(length(epi_rank) > 0){
          rkd_epiunits$candidates[rkd_epiunits$source == sourcenodes[a]][1:length(epi_rank)] <- as.integer(names(head_of(networks_2monthsep[[i]], edgz)[epi_rank])) #puts them in order, the rest will be NAs
        }
        edgz <- epi_rank <- NA
      }
      
      #figure out what the largest out-degree of the source node is, giving us the number of ranks to look through
      #mostedges <- max(degree(networks_2monthsep[[i]],v=as.character(sourcenodes),mode = "out"))
      
      #need to keep track of the total nodes being searched overall
      if(j == 1){
        totalnodes <- rkd_epiunits$candidates[rkd_epiunits$rank == 1]
        totalnodes <- unique(totalnodes)
        totalnotes <- totalnodes[!is.na(totalnodes)]
      }
      if(j > 1){
        totalnodes <- append(totalnodes,rkd_epiunits$candidates[rkd_epiunits$rank == 1][which(!(rkd_epiunits$candidates[rkd_epiunits$rank == 1] %in% totalnodes))])
        totalnodes <- unique(totalnodes)
        totalnotes <- totalnodes[!is.na(totalnodes)]
      }
      for(m in 1:mostedges){
        if(m > 1){ #keep track of how many of the epiunits have been searched
          totalnodes <- append(totalnodes,rkd_epiunits$candidates[rkd_epiunits$rank == m][which(!(rkd_epiunits$candidates[rkd_epiunits$rank == m] %in% totalnodes))])
          totalnodes <- unique(totalnodes)
          totalnodes <- totalnodes[!is.na(totalnodes)]
        }
        #for all of the source nodes, look at their mth ranked-edge-epiunit and record the number of epiunits looked at and whether any of them are t+1 outbreaks
        #check whether the number of nodes looked at is >= x%
        #fill in all of the information into wctracing_outneighbs and wcontacttracing_tplus1epiunits
        #once >= 35%, break out of this loop
        
        for(l in 1:length(pct_smol)){
          if(l > counter){ #go through only the radii that haven't been filled in yet
            if((length(totalnodes)) >= (numtotalepiunits*pct_smol[l])){ #if see greater than or = the number of epiunits at that percentage, fill in the information
              
              wctracing_outneighbs$pctntwk[wctracing_outneighbs$twomthntwk == i & wctracing_outneighbs$pctaway == pct_smol[l]] <- length(totalnodes)/numtotalepiunits
              
              wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$twomthntwk == i & wctracing_outneighbs$pctaway == pct_smol[l]]  <- 
                length(totalnodes[totalnodes %in% tplus1outbreaks])/length(outbreakepiunits[[(i+1)]])
              
              wctracing_outneighbs$edgerank[wctracing_outneighbs$twomthntwk == i & wctracing_outneighbs$pctaway == pct_smol[l]] <- m
              
              wctracing_outneighbs$order[wctracing_outneighbs$twomthntwk == i & wctracing_outneighbs$pctaway == pct_smol[l]] <- j
              
              wcontacttracing_tplus1epiunits[[((i-1)*(length(pct_smol)) + l)]] <- totalnodes[totalnodes %in% tplus1outbreaks]
              
              counter <- l #if go into this if statement, mean you're > goal % age for pct[l] so don't want to go in there again
              print(paste("counter =", counter, "order = ", j))
              if(counter == length(pct_smol)){
                sourcenodes <- mostedges <- rkd_epiunits <- NA
                break
              }
            }
          }
        }
        #break out of the m loop
        if(counter == length(pct_smol)){
          sourcenodes <- mostedges <- rkd_epiunits <- NA
          break
        } #don't need to keep going through the steps, this only breaks you out of the forloop it's in..not the whole nested set
        
      }
      #once have searched >= 35%, break out of the j loop also
      if(counter == length(pct_smol)){
        sourcenodes <- mostedges <- rkd_epiunits <- NA
        break
      } #don't need to keep going through the steps, this only breaks you out of the forloop it's in..not the whole nested set
      
      #figure out what the sourcenodes are for (j+1)th order
      newsources <- unique(rkd_epiunits$candidates)
      newsources <- newsources[!is.na(newsources)]
      
      sourcenodes <- mostedges <- rkd_epiunits <- NA
      
    }
    counter <- NULL
    mostedges <- sourcenodes <- rkd_epiunits <- tplus1outbreaks <- NULL
  }
  end_time <- Sys.time()
  end_time - start_time 
  
  #saveRDS(wctracing_outneighbs, file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/SerotypeSpecific/Serotype",serotype,"/2months_weightedcontacttracing_outneighbourhoodmetrics_11.26.2023.rds"))
  #saveRDS(wcontacttracing_tplus1epiunits, file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/SerotypeSpecific/Serotype",serotype,"/2months_weightedcontacttracingtplus1epiunits_11.26.2023.rds"))
  
  
  outbreakepiunits <- tplus1stats_full_fxd <- t_outbreaks <- serotype <- wctracing_outneighbs <- wcontacttracing_tplus1epiunits <- NULL
}

####MAKING TABLE S3, S4####
#need outbreakepiunits_A/_O, t_outbreaks_a/_o (because only looking at the t outbreaks that act as sources in the network), tplus1_stats_a/_o (because only going to find outbreaks that have network destination data), outbreakepiunits

numserotypes <- 2
numtotalepiunits <- 54096
p_perc = c(0.05,0.1,0.15,0.20,0.25,0.3,0.35)

for(n in 1:numserotypes){
  
  if(n == 1){
    outbreakepiunits <- outbreakepiunits_A
    tplus1stats_full_fxd <- tplus1_stats_a
    t_outbreaks <- t_outbreaks_a
    serotype <- "A"
  }
  
  if(n == 2){
    outbreakepiunits <- outbreakepiunits_O
    tplus1stats_full_fxd <- tplus1_stats_o
    t_outbreaks <- t_outbreaks_o
    serotype <- "O"
  }
  
  #pperc_outbreaks_full <- readRDS(file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/SerotypeSpecific/Serotype",serotype,"/2monthnetworks_degreeonlyriskfactor_fullnetworkvalues_11.26.2023.rds")) 
  pperc_outbreaks_only2months <- pperc_outbreaks_full[pperc_outbreaks_full$twomonthpd != "alltime",]
  #euclideandist_varr_outneighbs <- readRDS(file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/SerotypeSpecific/Serotype",serotype,"/2months_euclidean_varr_distoutneighbourhoodmetrics_11.26.2023.rds")) 
  #wctracing_outneighbs <- readRDS(file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/SerotypeSpecific/Serotype",serotype,"/2months_weightedcontacttracing_outneighbourhoodmetrics_11.26.2023.rds"))
  
  for(k in 1:length(p_perc)){
    print(paste(p_perc[k],"% and serotype", serotype))
    print("risk factor (degree only) 2 months")
    a <- mean(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$P_perc == p_perc[k]], na.rm = T)
    print(a)
    print("risk factor (degree only) all time")
    b <- mean(pperc_outbreaks_full$Actual_peroutbreaks[pperc_outbreaks_full$P_perc == p_perc[k] & pperc_outbreaks_full$twomonthpd == "alltime"], na.rm = T)
    print(b)
    print("euclidean distance variable radius")
    c <- mean(euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$pctaway == p_perc[k]], na.rm = T)
    print(c)
    print("weighted contact tracing")
    d <- mean(wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$pctaway == p_perc[k]], na.rm = T)
    print(d)
    a <- b <- c <- d <- NULL
  }
  
}
