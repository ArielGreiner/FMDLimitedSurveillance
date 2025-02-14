#######Calculate spatial proximity, network connectivity (2 months), network connectivity (all time) and network proximity at 40-100%

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

####LOAD IN USEFUL FUNCTIONS####

#functions written by José L. Herrera Diestra

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

p_perc_big <- c(0.05,0.1,0.15,0.20,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)


####SPATIAL PROXIMITY METHOD####

#need to load in: outbreakData, outbreakepiunits, locationsmasterlist
numtotalepiunits <- 54096 #dim(locations)[1] 

#run a for loop where, for every month, the radius slowly increases until it exceeds x%*(totalepiunits) and then it records data from that and continues until it hits 100% of all epiunits
steps <- 700 #need to change this to 700 because twomonthntwk 22 only got up to 80% if only let it get up to 500, everything else seems fine
pct <- c(0.05,0.1,0.15,0.20,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)

euclideandist_varr_outneighbs <- data.frame(twomthntwk = rep(seq(1,num_2_ntwks,1), each = length(pct)), pctaway = rep(pct, num_2_ntwks), pctntwk = NA, avgsize_peroutbk = NA, pcttplus1outbks = NA, radii = NA)

proximity_varr_tplus1epiunits <- list()

for(i in 1:num_2_ntwks){
  print(paste("i = ", i))
  t_outbreakdata <- locationsmasterlist[(locationsmasterlist$vill_ID %in% outbreakepiunits[[i]]),]
  counter <- 0
  for(j in 1:steps){
    stepsize <- j/50 #chosen by trial and error, aiming for a stepsize that was small enough that the resulting surveillance effort %age was close to what we were aiming for, see Fig. S1
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

#save euclideandist_varr_outneighbs and proximity_varr_tplus1epiunits, see below for how i saved them
#saveRDS(euclideandist_varr_outneighbs, file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/Beyond35/2months_euclidean_varr_distoutneighbourhoodmetrics_1.2025.rds")) 
#saveRDS(proximity_varr_tplus1epiunits, file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/Beyond35/2months_proximityvarr_tplus1epiunits_1.2025.rds")) 


####NETWORK CONNECTIVITY METHODS (2 MONTHS + ALL TIME)####

#need to load in networkmetrics_2monthsep, outbreakepiunits, epiunit_networkmetrics, outbreakData
#How well does the 2-month network predict the 2nd-month outbreaks? (t,t+1) -> (t+1) version

#divide by the number of outbreaks observed
p_perc = c(0.05,0.1,0.15,0.20,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)
metrix <- c("Degree")
numtotalepiunits <- 54096

pperc_outbreaks_full <- data.frame(metric = rep(rep(metrix, each = length(p_perc)),(num_2_ntwks+1)),P_perc = rep(rep(p_perc,length(metrix)), (num_2_ntwks+1)), twomonthpd = rep(c(seq(1,num_2_ntwks,1), "alltime"), each = (length(p_perc)*length(metrix))), Actual_peroutbreaks = NA)


for(j in 1:num_2_ntwks){
  topDegreeepiunits <- list()
  networkmetrics_2monthsep_indiv <- networkmetrics_2monthsep[[j]] 
  #num_epiunits <- gorder(networks_2monthsep[[j]]) #removed 9.28.2023 so i could put these values on the same scale as the euclidean distance and contact tracing full network methods
  onemonthoutbreaks <- outbreakepiunits[[(j+1)]]
  
  for(i in 1:length(p_perc)){
    #which are the top P% of nodes for each metric?
    #note that if multiple nodes have the same degree, their order is determined by the order that the network creation function (networkCreationDW) ordered them in (which is itself determined by igraph)
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
outbreakepiunits_nontemp <- unique(outbreakData$vill_ID)

#Actual = any epiunits that had an outbreak at ANY point
#divide by the number of outbreaks observed
metrix <- c("Degree")
pperc_outbreaks <- data.frame(metric = rep(metrix, each = length(p_perc)),P_perc = rep(p_perc,length(metrix)), Actual_peroutbreaks = NA)

for(i in 1:length(p_perc)){
  pperc_outbreaks$Actual_peroutbreaks[pperc_outbreaks$metric == "Degree" & pperc_outbreaks$P_perc == p_perc[i]] <- length(which(topDegreeepiunits[[i]] %in% outbreakepiunits_nontemp))/length(outbreakepiunits_nontemp) 
  
}

#add it to the temporal networks dataframes
pperc_outbreaks_full$Actual_peroutbreaks[pperc_outbreaks_full$twomonthpd == "alltime"] <- pperc_outbreaks$Actual_peroutbreaks

#save pperc_outbreaks_full, see below for how I saved them
#saveRDS(pperc_outbreaks_full, file = "FMDLimitedSurveillance/TemporalNetworks_7.2023/Beyond35/2monthnetworks_degreeonlyriskfactor_fullnetworkvalues_1.2025.rds") 

####NETWORK PROXIMITY METHOD####

#need networks_2monthsep, t_outbreaks (because only looking at the t outbreaks that act as sources in the network), tplus1stats_full_fxd (because only going to find outbreaks that have network destination data), outbreakepiunits

numtotalepiunits <- 54096 
pct_smol <- c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)
maxorder = 20

wctracing_outneighbs <- data.frame(twomthntwk = rep(seq(1,num_2_ntwks,1), each = length(pct_smol)), pctaway = rep(pct_smol, num_2_ntwks), pctntwk = NA, pcttplus1outbks = NA, edgerank = NA, order = NA)
wcontacttracing_tplus1epiunits <- list()

start_time <- Sys.time()
for(i in 1:num_2_ntwks){
  print(paste("i = ", i))
  tplus1outbreaks <- tplus1stats_full_fxd$epiunit[tplus1stats_full_fxd$inf_status == 1 & tplus1stats_full_fxd$tplus1 == (i+1)]
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
      srchdnodes <- append(srchdnodes, sourcenodes)
      newsources <- NA
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
      epi_rank <- order(edgz$weight, decreasing = TRUE) #for ties, it goes in the order that they were listed (determined by igraph)
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

#need to save wctracing_outneighbs and wcontacttracing_tplus1epiunits, showing how i did it below
#saveRDS(wctracing_outneighbs, file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/Beyond35/2months_weightedcontacttracing_outneighbourhoodmetrics_1.2025.rds"))
#saveRDS(wcontacttracing_tplus1epiunits, file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/Beyond35/2months_weightedcontacttracingtplus1epiunits_1.2025.rds"))


####PLOTTING EVERYTHING####
t_outbreaks <- readRDS(file = paste0("FMDLimitedSurveillance/Data/t_outbreaks.rds"))
tplus1stats_full_fxd <- readRDS(file = paste0("FMDLimitedSurveillance/Data/tplus1stats_full_fxd.rds"))


#data to be loaded in:
numtotalepiunits <- 54096
p_perc = c(0.05,0.1,0.15,0.20,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)

pperc_outbreaks_full <- readRDS(file = "FMDLimitedSurveillance/TemporalNetworks_7.2023/Beyond35/2monthnetworks_degreeonlyriskfactor_fullnetworkvalues_1.2025.rds") 
pperc_outbreaks_only2months <- pperc_outbreaks_full[pperc_outbreaks_full$twomonthpd != "alltime",]

euclideandist_varr_outneighbs <- readRDS(file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/Beyond35/2months_euclidean_varr_distoutneighbourhoodmetrics_1.2025.rds")) 
wctracing_outneighbs <- readRDS(file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/Beyond35/2months_weightedcontacttracing_outneighbourhoodmetrics_1.2025.rds"))

network_tplus1outbreaks <- rep(NA,num_2_ntwks)
tplus1outbreaks <- rep(NA, num_2_ntwks)
#pperc_outbreaks_avgmetrics <- data.frame(P_perc = rep(p_perc, num_2_ntwks), twomonthpd = rep(seq(1,num_2_ntwks,1), each = length(p_perc)), Actual_peroutbreaks = NA) #dont need this anymore bc just one metric

for(i in 1:num_2_ntwks){
  network_tplus1outbreaks[i] <- length(tplus1stats_full_fxd$epiunit[tplus1stats_full_fxd$inf_status == 1 & tplus1stats_full_fxd$tplus1 == (i+1)])
  tplus1outbreaks[i] <- length(outbreakepiunits[[(i+1)]])
  #for(j in 1:length(p_perc)){
  # pperc_outbreaks_avgmetrics$Actual_peroutbreaks[pperc_outbreaks_avgmetrics$P_perc == p_perc[j] & pperc_outbreaks_avgmetrics$twomonthpd == i] <- pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$P_perc == p_perc[j] & pperc_outbreaks_only2months$twomonthpd == i]
  #}
}

#random method - re-run the code
#10000 replicates

set.seed(2)
numreps <- 10000
randomsamps <- list()
random_numoutbks_actual_peroutbreaks <- matrix(data = NA, nrow = numreps, ncol = length(p_perc))

random_method <- data.frame(P_perc = rep(p_perc), mean_peroutbreaks = NA, var_peroutbreaks = NA, twentyfive_peroutbreaks = NA, seventyfive_peroutbreaks = NA)
alloutbreaks <- unique(unlist(outbreakepiunits))

for(j in 1:numreps){
  randomsamps[[j]] <- sample(seq(1,numtotalepiunits), numtotalepiunits*max(p_perc)) #generate a sample that's the largest size needed
  for(i in 1:length(p_perc)){
    random_numoutbks_actual_peroutbreaks[j,i] <- length(which(alloutbreaks %in% randomsamps[[j]][1:(numtotalepiunits*p_perc[i])]))/length(alloutbreaks)
  }
}
for(i in 1:length(p_perc)){
  random_method$mean_peroutbreaks[random_method$P_perc == p_perc[i]] <- mean(random_numoutbks_actual_peroutbreaks[,i])
  random_method$var_peroutbreaks[random_method$P_perc == p_perc[i]] <- var(random_numoutbks_actual_peroutbreaks[,i])
  random_method$twentyfive_peroutbreaks[random_method$P_perc == p_perc[i]] <- quantile(random_numoutbks_actual_peroutbreaks[,i], 0.25) 
  random_method$seventyfive_peroutbreaks[random_method$P_perc == p_perc[i]] <- quantile(random_numoutbks_actual_peroutbreaks[,i], 0.75) 
}


#take the averages of all of the different methods and then put it into one dataframe
avgresults <- data.frame(P_perc = rep(p_perc), Random = NA, Proximity = NA, RiskFactors = NA, ContactTracing = NA)

avgresults$Random <- random_method$mean_peroutbreaks

#changing it so that when there is an NA in a month, it shows up as an NA so can plot 5-100% for all of them
#but need to remove month 66 because it's wrong
for(i in 1:length(p_perc)){
  avgresults$Proximity[avgresults$P_perc == p_perc[i]] <- mean(euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$pctaway == p_perc[i] & euclideandist_varr_outneighbs$twomthntwk < 66]) #, na.rm = T
  avgresults$RiskFactors[avgresults$P_perc == p_perc[i]] <- mean(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$P_perc == p_perc[i] & pperc_outbreaks_only2months$twomonthpd < 66]) #, na.rm = T
  avgresults$ContactTracing[avgresults$P_perc == p_perc[i]] <- mean(wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$pctaway == p_perc[i] & wctracing_outneighbs$twomthntwk < 66]) #, na.rm = T
}

#avgresults_abr <- avgresults[avgresults$P_perc %in% c(0.05,0.15,0.30),] #not using this for this set of results

#6.11.2024: take the averages of all of the different methods and then put it into one dataframe + 25%, 75% quantiles
avgresults_exp <- data.frame(P_perc = rep(p_perc*100, each = 4), SearchMethod = c("Random", "Proximity", "RiskFactors", "ContactTracing"), Mean = NA, TwentyFive = NA, SeventyFive = NA)

avgresults_exp$Mean[avgresults_exp$SearchMethod == "Random"] <- random_method$mean_peroutbreaks
avgresults_exp$TwentyFive[avgresults_exp$SearchMethod == "Random"] <- random_method$twentyfive_peroutbreaks
avgresults_exp$SeventyFive[avgresults_exp$SearchMethod == "Random"] <- random_method$seventyfive_peroutbreaks

for(i in 1:length(p_perc)){ #had to run this one by one because of the error messages, just because when it hit NAs it broke out of the loop
  avgresults_exp$Mean[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "Proximity"] <- mean(euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$pctaway == p_perc[i] & euclideandist_varr_outneighbs$twomthntwk < 66]) #, na.rm = T
  
  avgresults_exp$TwentyFive[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "Proximity"] <- quantile(euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$pctaway == p_perc[i] & euclideandist_varr_outneighbs$twomthntwk < 66], 0.25) #, na.rm = T
  
  avgresults_exp$SeventyFive[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "Proximity"] <- quantile(euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$pctaway == p_perc[i] & euclideandist_varr_outneighbs$twomthntwk < 66], 0.75) #, na.rm = T
  
  if(p_perc[i] < 0.65){
    avgresults_exp$Mean[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "RiskFactors"] <- mean(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$P_perc == p_perc[i] & pperc_outbreaks_only2months$twomonthpd < 66]) #, na.rm = T
    
    avgresults_exp$TwentyFive[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "RiskFactors"] <- quantile(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$P_perc == p_perc[i] & pperc_outbreaks_only2months$twomonthpd < 66], 0.25) #, na.rm = T
    
    avgresults_exp$SeventyFive[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "RiskFactors"] <- quantile(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$P_perc == p_perc[i] & pperc_outbreaks_only2months$twomonthpd < 66], 0.75) #, na.rm = T
  }
  if(p_perc[i] < 0.35){
    avgresults_exp$Mean[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "ContactTracing"] <- mean(wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$pctaway == p_perc[i] & wctracing_outneighbs$twomthntwk < 66]) #, na.rm = T
    
    avgresults_exp$TwentyFive[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "ContactTracing"] <- quantile(wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$pctaway == p_perc[i] & wctracing_outneighbs$twomthntwk < 66], 0.25) #, na.rm = T
    
    avgresults_exp$SeventyFive[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "ContactTracing"] <- quantile(wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$pctaway == p_perc[i] & wctracing_outneighbs$twomthntwk < 66], 0.75) #, na.rm = T
  }
}


#6.11.2024: Plot avgresults_exp (error bars!)
#7.7.2024 - multiply by 100
#1.31.2025: stopped network connectivity at 60% because some of the months aren't able to search beyond 63%; also adjusted the alpha and size of the points/lines
plot <- ggplot(avgresults_exp, aes(x = P_perc, y = Mean*100, ymin = TwentyFive*100, ymax = SeventyFive*100, fill = SearchMethod, col = SearchMethod)) +
  geom_ribbon(alpha = 0.3, colour = NA)+ #0.1
  geom_line(lwd = 1)+ #1.5
  geom_point(aes(shape=SearchMethod), size = 2.5)+ #3
  scale_shape_manual(values=c(15,17,18,16))+
  #geom_point()+
  xlab("Surveillance Effort")+
  ylab("% Outbreaks Detected")+
  scale_y_continuous(limits = c(0, 100))+
  scale_x_discrete(limits = p_perc*100)+
  scale_fill_manual(values=c("#beaed4", "#fdc086", "#7fc97f", "#386cb0"), breaks = c("Proximity", "RiskFactors", "ContactTracing", "Random"), labels = c("Proximity", "RiskFactors", "ContactTracing", "Random"), name = "Surveillance Method")+
  scale_colour_manual(values=c("#beaed4", "#fdc086", "#7fc97f", "#386cb0"), breaks = c("Proximity", "RiskFactors", "ContactTracing", "Random"), labels = c("Proximity", "RiskFactors", "ContactTracing", "Random"), name = "Surveillance Method")+
  theme_bw()
ggsave(plot, filename = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/Beyond35/1.2025_RevisionsGraphs/AvgMethods_ErrorBars_1.31.2025.png"), bg = "transparent", height = 10, width = 10) #1.16.2025

