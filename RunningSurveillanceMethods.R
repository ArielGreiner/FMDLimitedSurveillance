####CODE TO ASSESS THE SURVEILLANCE METHODS AND GENERATE THE FIGURES IN THE MANUSCRIPT AND SUPPLEMENTARY MATERIAL

#might not need all of these...
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

#functions written by José L. Herrera-Diestra

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
#need to load in: outbreakData, outbreakepiunits, locationsmasterlist
numtotalepiunits <- 54096 #dim(locations)[1] 

#run a for loop where, for every month, the radius slowly increases until it exceeds x%*(totalepiunits) and then it records data from that and continues until it hits >35% of all epiunits
steps <- 500
pct <- c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.5,0.75,1)

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
#saveRDS(euclideandist_varr_outneighbs, file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2months_euclidean_varr_distoutneighbourhoodmetrics_10.26.2023.rds")) 
#saveRDS(proximity_varr_tplus1epiunits, file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2months_proximityvarr_tplus1epiunits_10.26.2023.rds")) 



####NETWORK CONNECTIVITY METHOD####
#need to load in networkmetrics_2monthsep, outbreakepiunits, epiunit_networkmetrics, outbreakData
#How well does the 2-month network predict the 2nd-month outbreaks? (t,t+1) -> (t+1) version

#divide by the number of outbreaks observed
p_perc = c(0.05,0.1,0.15,0.20,0.25,0.3,0.35)
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
#saveRDS(pperc_outbreaks_full, file = "FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthnetworks_degreeonlyriskfactor_fullnetworkvalues_11.9.2023.rds") 

#Fig S3b: how does it change over time? remove the "alltime" rows first
pperc_outbreaks_only2months <- pperc_outbreaks_full[pperc_outbreaks_full$twomonthpd != "alltime",]

plot <- ggplot(data =pperc_outbreaks_only2months[pperc_outbreaks_only2months$P_perc==0.1,], aes(x=as.numeric(twomonthpd), y=Actual_peroutbreaks, group = metric))+
  geom_line(aes(color = metric))+
  ggtitle("Change in Metric Over Time")+
  #scale_x_discrete(limits = p_perc)+
  xlab("2 Month Period")+
  ylab("% Outbreaks Happening in Top 10% of Nodes")+
  ylim(c(0,1)) #added 7.2.2024
#scale_color_viridis(discrete = TRUE)
#scale_color_manual(values = c(sixtyfiveviridiscols, "red"))

#Fig S3a (first row): if look at the top x% of nodes, how many t+1 outbreaks do we find? looking at the degree metric and across networks (e.g. months)

pperc_outbreaks_avg <- data.frame(P_perc = p_perc, avg_tplus1outbreaks = NA, avg_alloutbreaks = NA, lower_tplus1outbreaks = NA, lower_alloutbreaks = NA, higher_tplus1outbreaks = NA, higher_alloutbreaks = NA)

#removing month (66,67) because no outbreaks observed in month 67
pperc_outbreaks_only2months <- pperc_outbreaks_only2months[pperc_outbreaks_only2months$twomonthpd < 66,]

for(i in 1:length(p_perc)){
  pperc_outbreaks_avg$avg_tplus1outbreaks[pperc_outbreaks_avg$P_perc == p_perc[i]] <- mean(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$P_perc == p_perc[i]])
  
  pperc_outbreaks_avg$lower_tplus1outbreaks[pperc_outbreaks_avg$P_perc == p_perc[i]] <- range(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$P_perc == p_perc[i]])[1]
  
  pperc_outbreaks_avg$higher_tplus1outbreaks[pperc_outbreaks_avg$P_perc == p_perc[i]] <- range(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$P_perc == p_perc[i]])[2]
  
  pperc_outbreaks_avg$avg_alloutbreaks[pperc_outbreaks_avg$P_perc == p_perc[i]] <- mean(pperc_outbreaks_full$Actual_peroutbreaks[pperc_outbreaks_full$P_perc == p_perc[i] & pperc_outbreaks_full$twomonthpd == "alltime"])
  
  pperc_outbreaks_avg$lower_alloutbreaks[pperc_outbreaks_avg$P_perc == p_perc[i]] <- range(pperc_outbreaks_full$Actual_peroutbreaks[pperc_outbreaks_full$P_perc == p_perc[i] & pperc_outbreaks_full$twomonthpd == "alltime"])[1]
  
  pperc_outbreaks_avg$higher_alloutbreaks[pperc_outbreaks_avg$P_perc == p_perc[i]] <- range(pperc_outbreaks_full$Actual_peroutbreaks[pperc_outbreaks_full$P_perc == p_perc[i] & pperc_outbreaks_full$twomonthpd == "alltime"])[2]
}

#track which epiunits were searched for fig S8
riskfactor_degreeonly_tplus1epiunits <- list()

#divide by the number of outbreaks observed
p_perc = c(0.05,0.1,0.15,0.20,0.25,0.3,0.35)
numtotalepiunits <- 54096

for(j in 1:num_2_ntwks){
  networkmetrics_2monthsep_indiv <- networkmetrics_2monthsep[[j]] 
  onemonthoutbreaks <- outbreakepiunits[[(j+1)]]
  
  for(i in 1:length(p_perc)){
    #which are the top P% of nodes for each metric?
    itera <- ((j-1)*(length(p_perc)) + i)
    riskfactor_degreeonly_tplus1epiunits[[itera]] <- networkmetrics_2monthsep_indiv$Node[order(networkmetrics_2monthsep_indiv$Degree, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)][which(networkmetrics_2monthsep_indiv$Node[order(networkmetrics_2monthsep_indiv$Degree, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)] %in% onemonthoutbreaks)]
    itera <- NULL
  }
}
#need to save riskfactor_degreeonly_tplus1epiunits, how I did it is shown below
#saveRDS(riskfactor_degreeonly_tplus1epiunits, file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2months_riskfactortplus1epiunits_10.14.2023.rds")) 



####NETWORK PROXIMITY METHOD####
#need networks_2monthsep, t_outbreaks (because only looking at the t outbreaks that act as sources in the network), tplus1stats_full_fxd (because only going to find outbreaks that have network destination data), outbreakepiunits

numtotalepiunits <- 54096 
pct_smol <- c(0.05,0.1,0.15,0.2,0.25,0.3,0.35)
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
#saveRDS(wctracing_outneighbs, file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2months_weightedcontacttracing_outneighbourhoodmetrics_11.5.2023.rds"))
#saveRDS(wcontacttracing_tplus1epiunits, file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2months_weightedcontacttracingtplus1epiunits_11.5.2023.rds"))



####COMPARING THE THREE METHODS####

#data to be loaded in: outbreakepiunits, networks_2monthsep, tplus1stats_full_fxd

#and also, from network connectivity method: pperc_outbreaks_full, pperc_outbreaks_only2months 
#pperc_outbreaks_full <- readRDS(file = "FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthnetworks_degreeonlyriskfactor_fullnetworkvalues_11.9.2023.rds") 
#pperc_outbreaks_only2months <- pperc_outbreaks_full[pperc_outbreaks_full$twomonthpd != "alltime",]

#and also, from the spatial proximity method: euclideandist_varr_outneighbs
#euclideandist_varr_outneighbs <- readRDS(file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2months_euclidean_varr_distoutneighbourhoodmetrics_10.26.2023.rds"))

#and lastly, from the network proximity method: wctracing_outneighbs
#wctracing_outneighbs <- readRDS(file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2months_weightedcontacttracing_outneighbourhoodmetrics_11.5.2023.rds"))

numtotalepiunits <- 54096
#also p_perc 
p_perc = c(0.05,0.1,0.15,0.20,0.25,0.3,0.35)


network_tplus1outbreaks <- rep(NA,num_2_ntwks)
tplus1outbreaks <- rep(NA, num_2_ntwks)

for(i in 1:num_2_ntwks){
  network_tplus1outbreaks[i] <- length(tplus1stats_full_fxd$epiunit[tplus1stats_full_fxd$inf_status == 1 & tplus1stats_full_fxd$tplus1 == (i+1)])
  tplus1outbreaks[i] <- length(outbreakepiunits[[(i+1)]])
}


for(k in 1:length(p_perc)){ 
  
  #Fig S2
  #png(paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthplots/ThreeMethodComparison_11.2023/searchsize_Comparing3Methods_",p_perc[k]*100,"pct_all2monthnetworks_11.9.2023.png"))
  plot(x = seq(1,num_2_ntwks,1), y = (wctracing_outneighbs$pctntwk[wctracing_outneighbs$pctaway == p_perc[k]]*100), xlab = "Month", ylab = "% Network Searched", pch = 20, col = "#7fc97f", main = paste("Searching ~",p_perc[k]*100,"% of Network"), ylim = c(0,60))
  lines(x = seq(1,num_2_ntwks,1), y = (wctracing_outneighbs$pctntwk[wctracing_outneighbs$pctaway == p_perc[k]]*100), lty = 2, col = "#7fc97f")
  points(x = seq(1,num_2_ntwks,1), y = (euclideandist_varr_outneighbs$pctntwk[euclideandist_varr_outneighbs$pctaway == p_perc[k]]*100), pch=20, col = "#beaed4")
  lines(x = seq(1,num_2_ntwks,1), y = (euclideandist_varr_outneighbs$pctntwk[euclideandist_varr_outneighbs$pctaway == p_perc[k]]*100), lty=2, col = "#beaed4")
  points(x = seq(1,num_2_ntwks,1), y = rep(p_perc[k]*100,num_2_ntwks), pch=20, col = "#fdc086")
  lines(x = seq(1,num_2_ntwks,1), y = rep(p_perc[k]*100,num_2_ntwks), lty=2, col = "#fdc086")
  legend(1, 60, legend=c("Proximity (Var r)", "Risk Factors","Weighted Contact Tracing"),col=c("#beaed4","#fdc086","#7fc97f"), pch = rep(20,3), cex=0.8)
  #dev.off()
  
  
  #Make some plots showing which method is the best in each month, characterized a few different ways, for 5/10/15/20/25/30/35%
  #first step - make a dataframe that contains the number of t+1 outbreaks found in each different way, for 5/10/15/20/25/30/35%
  #then add a column giving the # of t+1 outbreaks in said month (i = i+1 in this setting)
  #then add the size of network searched each different way, for 5/10/15/20/25/30/35%
  #then add a column saying which method found the most t+1 outbreaks and another column saying which had the highest 't+1 outbreaks found/amount of network searched' ratio
  summary_threemethods <- data.frame(tplus1month = seq(1,num_2_ntwks,1), tplus1outbreaks = tplus1outbreaks, tplus1_proximity = (euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$pctaway == p_perc[k]]*tplus1outbreaks), tplus1_riskfactors = (pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$P_perc == p_perc[k]]*tplus1outbreaks), tplus1_ctracing = (wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$pctaway == p_perc[k]]*tplus1outbreaks), 
                                     pctntwk_proximity = (euclideandist_varr_outneighbs$pctntwk[euclideandist_varr_outneighbs$pctaway == p_perc[k]]*numtotalepiunits), pctntwk_riskfactors = rep((p_perc[k]*numtotalepiunits),num_2_ntwks), pctntwk_ctracing = (wctracing_outneighbs$pctntwk[wctracing_outneighbs$pctaway == p_perc[k]])*numtotalepiunits, best_numoutbreaks_id = NA, best_numoutbreaks_val = NA, best_ratio_id = NA, best_ratio_val = NA) #removed 2.5.2024: tplus1outbreaks = tplus1outbreaks*p_perc[k]
  
  #following: https://stackoverflow.com/questions/43455862/add-new-column-with-name-of-max-column-in-data-frame and https://www.statology.org/r-max-across-columns/
  #10.31.2023 - changed so the ratio is a ratio of %-ages and not a ratio of #/# because this makes more sense
  summary_threemethods_dummy <- summary_threemethods[,c(3,4,5)]
  names(summary_threemethods_dummy) <- c("proximity", "riskfactors", "ctracing")
  summary_threemethods$best_numoutbreaks_id <- names(summary_threemethods_dummy)[max.col(summary_threemethods_dummy)]
  summary_threemethods$best_numoutbreaks_val <- pmax(summary_threemethods_dummy$proximity, summary_threemethods_dummy$riskfactors, summary_threemethods_dummy$ctracing)
  
  summary_threemethods_dummy <- data.frame(riskfactors = pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$P_perc == p_perc[k]]/(p_perc[k]), 
                                           proximity = euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$pctaway == p_perc[k]]/euclideandist_varr_outneighbs$pctntwk[euclideandist_varr_outneighbs$pctaway == p_perc[k]], 
                                           ctracing = (wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$pctaway == p_perc[k]])/((wctracing_outneighbs$pctntwk[wctracing_outneighbs$pctaway == p_perc[k]])))
  summary_threemethods$best_ratio_id <- names(summary_threemethods_dummy)[max.col(summary_threemethods_dummy)]
  summary_threemethods$best_ratio_val <- pmax(summary_threemethods_dummy$proximity, summary_threemethods_dummy$riskfactors, summary_threemethods_dummy$ctracing)
  
  summary_threemethods_descending <- summary_threemethods %>%
    arrange(desc(tplus1outbreaks))
  
 #fig 4 and s7
  summary_threemethods_descending_fixed <- summary_threemethods_descending[-66,]
  plot <- ggplot(summary_threemethods_descending_fixed, aes(x = reorder(tplus1month, -tplus1outbreaks), y = (best_numoutbreaks_val/tplus1outbreaks)*100), group = as.factor(best_numoutbreaks_id)) +
    geom_bar(aes(), fill = "lightgrey", stat = "identity")+
    geom_col(aes(y = (p_perc[k]*100), fill = as.factor(best_numoutbreaks_id)), col = "lightgrey")+ 
    #geom_col(aes(y = (p_perc[k]*100), fill = as.factor(best_numoutbreaks_id), col = as.factor(best_numoutbreaks_id)))+
    #geom_point(aes(y = 1, color = as.factor(best_numoutbreaks_id)))+
    #scale_fill_manual(values = c("lightgrey"))+
    scale_fill_manual(values = c("#7fc97f", "#beaed4","#fdc086"))+
    #scale_color_manual(values = c("#7fc97f", "#beaed4","#fdc086"))+
    geom_line(aes(x = reorder(tplus1month, -tplus1outbreaks), y = (p_perc[k]*100)), linewidth = 1, color = "black", group = 1)+
    xlab("Month (ordered by decreasing outbreak #)")+
    ylab("% t+1 Outbreaks")+
    scale_y_continuous(limits = c(0, 100))+
    ggtitle(paste("Method Finding Most t+1 Outbreaks, ",p_perc[k]*100,"% searched"))+
    theme_bw()
  #ggsave(plot, filename = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthplots/ThreeMethodComparison_11.2023/ThreeMethodSummary_highesttplus1outbreakpercentage_",p_perc[k]*100,"pct_greybars_3.14.2024.png"), bg = "transparent", height = 10, width = 10)
  

  #Fig 3, S6:Plot the # of outbreaks found by each method for each month, ordered by decreasing outbreak #
  plot <- ggplot(summary_threemethods_descending, aes(x = reorder(tplus1month, -tplus1outbreaks), y = tplus1_proximity)) +
    geom_bar(stat = "identity", fill = "#beaed4")+
    geom_line(aes(x = reorder(tplus1month, -tplus1outbreaks), y = tplus1outbreaks), linewidth = 1, color = "black", group = 1)+
    xlab("Month (ordered by decreasing outbreak #)")+
    ylab("% t+1 Outbreaks")+
    ggtitle(paste("% t+1 Outbreaks found by Proximity method, ",p_perc[k]*100,"% searched"))+
    theme_bw()
  #ggsave(plot, filename = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthplots/ThreeMethodComparison_11.2023/SingleMethodPlots/numproximitytplus1outbreaks_",p_perc[k]*100,"pct_1.4.2024.png"), bg = "transparent", height = 10, width = 10)
  
  plot <- ggplot(summary_threemethods_descending, aes(x = reorder(tplus1month, -tplus1outbreaks), y = tplus1_riskfactors)) +
    geom_bar(stat = "identity", fill = "#fdc086")+
    geom_line(aes(x = reorder(tplus1month, -tplus1outbreaks), y = tplus1outbreaks), linewidth = 1, color = "black", group = 1)+
    xlab("Month (ordered by decreasing outbreak #)")+
    ylab("% t+1 Outbreaks")+
    ggtitle(paste("% t+1 Outbreaks found by Risk Factor method, ",p_perc[k]*100,"% searched"))+
    theme_bw()
  #ggsave(plot, filename = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthplots/ThreeMethodComparison_11.2023/SingleMethodPlots/numriskfactorstplus1outbreaks_",p_perc[k]*100,"pct_1.4.2024.png"), bg = "transparent", height = 10, width = 10)
  
  plot <- ggplot(summary_threemethods_descending, aes(x = reorder(tplus1month, -tplus1outbreaks), y = tplus1_ctracing)) +
    geom_bar(stat = "identity", fill = "#7fc97f")+
    geom_line(aes(x = reorder(tplus1month, -tplus1outbreaks), y = tplus1outbreaks), linewidth = 1, color = "black", group = 1)+
    xlab("Month (ordered by decreasing outbreak #)")+
    ylab("% t+1 Outbreaks")+
    ggtitle(paste("% t+1 Outbreaks found by Contact Tracing method, ",p_perc[k]*100,"% searched"))+
    theme_bw()
  #ggsave(plot, filename = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthplots/ThreeMethodComparison_11.2023/SingleMethodPlots/numctracingtplus1outbreaks_",p_perc[k]*100,"pct_1.4.2024.png"), bg = "transparent", height = 10, width = 10)
  
  
  ##OUTBREAK OVERLAP
  #need to load in proximity_varr_tplus1epiunits, riskfactor_degreeonly_tplus1epiunits, wcontacttracing_tplus1epiunits
  #proximity_varr_tplus1epiunits <- readRDS(file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2months_proximityvarr_tplus1epiunits_10.26.2023.rds")) 
  #riskfactor_degreeonly_tplus1epiunits <- readRDS(file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2months_riskfactortplus1epiunits_10.14.2023.rds")) 
  #wcontacttracing_tplus1epiunits <- readRDS(file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2months_weightedcontacttracingtplus1epiunits_11.5.2023.rds"))
  
  p_perc = c(0.05,0.1,0.15,0.20,0.25,0.3,0.35)
  numstepsaway <- 20
  
  #now, compare these three (proximity_tplus1epiunits, wcontacttracing_tplus1epiunits, riskfactor_degreeonly_tplus1epiunits) with each other and with the risk factor ones
  
  #outbreakepiunits[[i+1]] are the t+1 outbreaks
  tplus1epiunit_overlap <- data.frame(tplus1month = rep(1,length(outbreakepiunits[[(1+1)]])), tplus1outbreaks = outbreakepiunits[[(1+1)]], proximity = 0, contacttracing = 0, riskfactor_degonly = 0, sum_degonly = 0)
  
  #for proximity, p_perc[k]*100%, so l = k in '((i-1)*(length(pct)) + k)' and i = tplus1month = 1
  tplus1epiunit_overlap$proximity[tplus1epiunit_overlap$tplus1outbreaks %in% proximity_varr_tplus1epiunits[[((1-1)*(length(pct)) + k)]]] <- 5
  #for weighted contact tracing, p_perc[k]*100% so l=k and j = tplus1month = 1 in "((j-1)*(length(pct_smol)) + l)"
  tplus1epiunit_overlap$contacttracing[tplus1epiunit_overlap$tplus1outbreaks %in% as.integer(wcontacttracing_tplus1epiunits[[((1-1)*(length(p_perc)) + k)]])] <- 100
  #for risk factor, p_perc[k]*100% is the kth element of p_perc so i = k and and j = tplus1month = 1 in '((j-1)*(length(p_perc)) + i)'
  tplus1epiunit_overlap$riskfactor_degonly[tplus1epiunit_overlap$tplus1outbreaks %in% riskfactor_degreeonly_tplus1epiunits[[((1-1)*(length(p_perc)) + k)]]] <- 10
  
  tplus1epiunit_overlap$sum_degonly <- tplus1epiunit_overlap$proximity + tplus1epiunit_overlap$contacttracing + tplus1epiunit_overlap$riskfactor_degonly
  
  for(j in 2:(num_2_ntwks-1)){
    dummy <- data.frame(tplus1month = rep(j,length(outbreakepiunits[[(j+1)]])), tplus1outbreaks = outbreakepiunits[[(j+1)]], proximity = 0, contacttracing = 0, riskfactor_degonly = 0, sum_degonly = 0)
    
    #for proximity, p_perc[k]*100%, so l = k in '((j-1)*(length(pct)) + k)' and j = tplus1month 
    dummy$proximity[dummy$tplus1outbreaks %in% proximity_varr_tplus1epiunits[[((j-1)*(length(pct)) + k)]]] <- 5
    #for weighted contact tracing, p_perc[k]*100% so l=k and j = tplus1month in "((j-1)*(length(pct_smol)) + l)"
    dummy$contacttracing[dummy$tplus1outbreaks %in% as.integer(wcontacttracing_tplus1epiunits[[((j-1)*(length(p_perc)) + k)]])] <- 100
    #for risk factor, p_perc[k]*100% is the kth element of p_perc so i = k and and j = tplus1month in '((j-1)*(length(p_perc)) + i)'
    dummy$riskfactor_degonly[dummy$tplus1outbreaks %in% riskfactor_degreeonly_tplus1epiunits[[((j-1)*(length(p_perc)) + k)]]] <- 10
    
    dummy$sum_degonly <- dummy$proximity + dummy$contacttracing + dummy$riskfactor_degonly
    
    tplus1epiunit_overlap <- rbind(tplus1epiunit_overlap, dummy)
  }
  
  
  #stacked bar chart ordered by # of outbreaks, each combo of methods gets a diff colour (rf-deg only)
  #need to make a new dataframe that summarizes, for each month, the number of epiunits that fall into each category
  #categories: "none", "prox", "rf", "rf+prox", "ct", "ct+prox", "ct+rf", "ct+rf+prox"
  methods <- c("none", "prox", "rf", "ct", "rf+prox", "ct+prox", "ct+rf", "ct+rf+prox")
  methods_fct <- factor(methods, levels = c("none", "prox", "rf", "ct", "rf+prox", "ct+prox", "ct+rf", "ct+rf+prox"))
  threemethodoverlap_permonth <- data.frame(tplus1month = rep(seq(1,(num_2_ntwks-1),1), each = length(methods)), method = rep(methods_fct), num_outbreaks_degonly = NA, num_outbreaks = NA)
  
  for(i in 1:(num_2_ntwks-1)){
    threemethodoverlap_permonth$num_outbreaks_degonly[threemethodoverlap_permonth$tplus1month == i & threemethodoverlap_permonth$method == "none"] <- length(tplus1epiunit_overlap$tplus1outbreaks[tplus1epiunit_overlap$tplus1month == i & tplus1epiunit_overlap$sum_degonly == 0])
    
    threemethodoverlap_permonth$num_outbreaks_degonly[threemethodoverlap_permonth$tplus1month == i & threemethodoverlap_permonth$method == "prox"] <- length(tplus1epiunit_overlap$tplus1outbreaks[tplus1epiunit_overlap$tplus1month == i & tplus1epiunit_overlap$sum_degonly == 5])
    
    threemethodoverlap_permonth$num_outbreaks_degonly[threemethodoverlap_permonth$tplus1month == i & threemethodoverlap_permonth$method == "rf"] <- length(tplus1epiunit_overlap$tplus1outbreaks[tplus1epiunit_overlap$tplus1month == i & tplus1epiunit_overlap$sum_degonly == 10])
    
    threemethodoverlap_permonth$num_outbreaks_degonly[threemethodoverlap_permonth$tplus1month == i & threemethodoverlap_permonth$method == "ct"] <- length(tplus1epiunit_overlap$tplus1outbreaks[tplus1epiunit_overlap$tplus1month == i & tplus1epiunit_overlap$sum_degonly == 100])
    
    threemethodoverlap_permonth$num_outbreaks_degonly[threemethodoverlap_permonth$tplus1month == i & threemethodoverlap_permonth$method == "rf+prox"] <- length(tplus1epiunit_overlap$tplus1outbreaks[tplus1epiunit_overlap$tplus1month == i & tplus1epiunit_overlap$sum_degonly == 15])
    
    threemethodoverlap_permonth$num_outbreaks_degonly[threemethodoverlap_permonth$tplus1month == i & threemethodoverlap_permonth$method == "ct+prox"] <- length(tplus1epiunit_overlap$tplus1outbreaks[tplus1epiunit_overlap$tplus1month == i & tplus1epiunit_overlap$sum_degonly == 105])
    
    threemethodoverlap_permonth$num_outbreaks_degonly[threemethodoverlap_permonth$tplus1month == i & threemethodoverlap_permonth$method == "ct+rf"] <- length(tplus1epiunit_overlap$tplus1outbreaks[tplus1epiunit_overlap$tplus1month == i & tplus1epiunit_overlap$sum_degonly == 110])
    
    threemethodoverlap_permonth$num_outbreaks_degonly[threemethodoverlap_permonth$tplus1month == i & threemethodoverlap_permonth$method == "ct+rf+prox"] <- length(tplus1epiunit_overlap$tplus1outbreaks[tplus1epiunit_overlap$tplus1month == i & tplus1epiunit_overlap$sum_degonly == 115])
  }
  

  #Fig S8: plot the number of outbreaks found by none of the methods in each month
  png(paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthplots/ThreeMethodComparison_11.2023/MethodEpiunitOverlap/",p_perc[k]*100,"pct/tplus1outbreaksfoundbynomethod_",p_perc[k]*100,"pct_month",i,"_degonly_7.2.2024.png"))
  plot(x = threemethodoverlap_permonth$tplus1month[threemethodoverlap_permonth$method == "none"], y = threemethodoverlap_permonth$num_outbreaks_degonly[threemethodoverlap_permonth$method == "none"]/tplus1outbreaks[-66], type = "l",xlab = "Month", ylab = "% of t+1 Outbreaks found by no method", main = paste(p_perc[k]*100,"% Searched, Degree Only"), ylim = c(0,1))
  points(x = threemethodoverlap_permonth$tplus1month[threemethodoverlap_permonth$method == "none"], y = threemethodoverlap_permonth$num_outbreaks_degonly[threemethodoverlap_permonth$method == "none"]/tplus1outbreaks[-66], pch=20)
  dev.off()
  
  threemethodoverlap_permonth <- summary_threemethods <- NA
}


####COMPARING THE THREE METHODS + RANDOM METHOD####
#load in the final results of all of the different methods: pperc_outbreaks_full, euclideandist_varr_outneighbs, wctracing_outneighbs
#pperc_outbreaks_full <- readRDS(file = "FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthnetworks_degreeonlyriskfactor_fullnetworkvalues_11.9.2023.rds") 
#euclideandist_varr_outneighbs <- readRDS(file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2months_euclidean_varr_distoutneighbourhoodmetrics_10.26.2023.rds"))
#wctracing_outneighbs <- readRDS(file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2months_weightedcontacttracing_outneighbourhoodmetrics_11.5.2023.rds"))

#also outbreakepiunits, networks_2monthsep, tplus1stats_full_fxd

#data to be loaded in:
numtotalepiunits <- 54096
#also p_perc 
p_perc = c(0.05,0.1,0.15,0.20,0.25,0.3,0.35)

pperc_outbreaks_fullonly <- pperc_outbreaks_full[pperc_outbreaks_full$twomonthpd == "alltime",]

pperc_outbreaks_only2months <- pperc_outbreaks_full[pperc_outbreaks_full$twomonthpd != "alltime",]

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


#take the averages of all of the different methods and then put it into one dataframe + 25%, 75% quantiles
avgresults_exp <- data.frame(P_perc = rep(p_perc*100, each = 4), SearchMethod = c("Random", "Proximity", "RiskFactors", "ContactTracing"), Mean = NA, TwentyFive = NA, SeventyFive = NA)

avgresults_exp$Mean[avgresults_exp$SearchMethod == "Random"] <- random_method$mean_peroutbreaks
avgresults_exp$TwentyFive[avgresults_exp$SearchMethod == "Random"] <- random_method$twentyfive_peroutbreaks
avgresults_exp$SeventyFive[avgresults_exp$SearchMethod == "Random"] <- random_method$seventyfive_peroutbreaks

for(i in 1:length(p_perc)){
  avgresults_exp$Mean[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "Proximity"] <- mean(euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$pctaway == p_perc[i]], na.rm = T)
  
  avgresults_exp$TwentyFive[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "Proximity"] <- quantile(euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$pctaway == p_perc[i]], 0.25, na.rm = T)
  
  avgresults_exp$SeventyFive[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "Proximity"] <- quantile(euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$pctaway == p_perc[i]], 0.75, na.rm = T)
  
  avgresults_exp$Mean[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "RiskFactors"] <- mean(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$P_perc == p_perc[i]], na.rm = T)
  
  avgresults_exp$TwentyFive[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "RiskFactors"] <- quantile(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$P_perc == p_perc[i]], 0.25, na.rm = T)
  
  avgresults_exp$SeventyFive[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "RiskFactors"] <- quantile(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$P_perc == p_perc[i]], 0.75, na.rm = T)
  
  avgresults_exp$Mean[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "ContactTracing"] <- mean(wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$pctaway == p_perc[i]], na.rm = T)
  
  avgresults_exp$TwentyFive[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "ContactTracing"] <- quantile(wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$pctaway == p_perc[i]], 0.25, na.rm = T)
  
  avgresults_exp$SeventyFive[avgresults_exp$P_perc == p_perc[i]*100 & avgresults_exp$SearchMethod == "ContactTracing"] <- quantile(wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$pctaway == p_perc[i]], 0.75, na.rm = T)
}


#Fig 2: Plot avgresults_exp (error bars!)
plot <- ggplot(avgresults_exp, aes(x = P_perc, y = Mean*100, ymin = TwentyFive*100, ymax = SeventyFive*100, fill = SearchMethod, col = SearchMethod)) +
  geom_ribbon(alpha = 0.3, colour = NA)+
  geom_line(lwd = 1)+
  geom_point(aes(shape=SearchMethod), size = 2.5)+
  scale_shape_manual(values=c(15,17,18,16))+
  xlab("Surveillance Effort")+
  ylab("% Outbreaks Detected")+
  scale_y_continuous(limits = c(0, 100))+
  scale_x_discrete(limits = p_perc*100)+
  scale_fill_manual(values=c("#beaed4", "#fdc086", "#7fc97f", "#386cb0"), breaks = c("Proximity", "RiskFactors", "ContactTracing", "Random"), labels = c("Proximity", "RiskFactors", "ContactTracing", "Random"), name = "Surveillance Method")+
  scale_colour_manual(values=c("#beaed4", "#fdc086", "#7fc97f", "#386cb0"), breaks = c("Proximity", "RiskFactors", "ContactTracing", "Random"), labels = c("Proximity", "RiskFactors", "ContactTracing", "Random"), name = "Surveillance Method")+
  theme_bw()
ggsave(plot, filename = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthplots/ThreeMethodComparison_11.2023/AvgMethods_ErrorBars_7.7.2024.png"), bg = "transparent", height = 10, width = 10)

#Fig S4: Plot avgresults_exp (error bars!) with all time network risk factor also
avgresults_exp_alltime <- data.frame(P_perc = rep(p_perc*100, each = 1), SearchMethod = "RiskFactors_AllTime", Mean = NA, TwentyFive = NA, SeventyFive = NA)
for(i in 1:length(p_perc)){ #just going to set 25% and 75% to the same value as 'Mean' so i can use the same format
  avgresults_exp_alltime$Mean[avgresults_exp_alltime$SearchMethod == "RiskFactors_AllTime" & avgresults_exp_alltime$P_perc == p_perc[i]*100] <- pperc_outbreaks_fullonly$Actual_peroutbreaks[pperc_outbreaks_fullonly$P_perc == p_perc[i]]
  avgresults_exp_alltime$TwentyFive[avgresults_exp_alltime$SearchMethod == "RiskFactors_AllTime" & avgresults_exp_alltime$P_perc == p_perc[i]*100] <- pperc_outbreaks_fullonly$Actual_peroutbreaks[pperc_outbreaks_fullonly$P_perc == p_perc[i]]
  avgresults_exp_alltime$SeventyFive[avgresults_exp_alltime$SearchMethod == "RiskFactors_AllTime" & avgresults_exp_alltime$P_perc == p_perc[i]*100] <- pperc_outbreaks_fullonly$Actual_peroutbreaks[pperc_outbreaks_fullonly$P_perc == p_perc[i]]
}

avgresults_exp_full <- rbind(avgresults_exp, avgresults_exp_alltime)

plot <- ggplot(avgresults_exp_full, aes(x = P_perc, y = Mean*100, ymin = TwentyFive*100, ymax = SeventyFive*100, fill = SearchMethod, col = SearchMethod)) +
  geom_ribbon(alpha = 0.3, colour = NA)+ #0.1
  geom_line(lwd = 1)+ #1.5
  geom_point(aes(shape=SearchMethod), size = 2.5)+ #3
  scale_shape_manual(values=c(15,17,18,20,16))+
  xlab("Surveillance Effort")+
  ylab("Average % Outbreaks Detected")+
  scale_y_continuous(limits = c(0, 100))+
  scale_x_discrete(limits = p_perc*100)+
  scale_fill_manual(values=c("#beaed4", "#fdc086", "#cc4c02", "#7fc97f", "#386cb0"), breaks = c("Proximity", "RiskFactors", "RiskFactors_AllTime", "ContactTracing", "Random"), labels = c("Proximity", "RiskFactors", "RiskFactors_AllTime", "ContactTracing", "Random"), name = "Surveillance Method")+
  scale_colour_manual(values=c("#beaed4", "#fdc086", "#cc4c02", "#7fc97f", "#386cb0"), breaks = c("Proximity", "RiskFactors", "RiskFactors_AllTime", "ContactTracing", "Random"), labels = c("Proximity", "RiskFactors", "RiskFactors_AllTime", "ContactTracing", "Random"), name = "Surveillance Method")+
  theme_bw()
ggsave(plot, filename = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthplots/ThreeMethodComparison_11.2023/AvgMethods_ErrorBars_Full_7.7.2024.png"), bg = "transparent", height = 10, width = 10)



####DETERMINING OTHER RELEVANT VALUES FOR THE MANUSCRIPT####
#data to be loaded in: pperc_outbreaks_full, euclideandist_varr_outneighbs, wctracing_outneighbs
#pperc_outbreaks_full <- readRDS(file = "FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthnetworks_degreeonlyriskfactor_fullnetworkvalues_11.9.2023.rds") 
#euclideandist_varr_outneighbs <- readRDS(file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2months_euclidean_varr_distoutneighbourhoodmetrics_10.26.2023.rds"))
#wctracing_outneighbs <- readRDS(file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2months_weightedcontacttracing_outneighbourhoodmetrics_11.5.2023.rds"))

numtotalepiunits <- 54096
#also p_perc 
p_perc = c(0.05,0.1,0.15,0.20,0.25,0.3,0.35)

pperc_outbreaks_only2months <- pperc_outbreaks_full[pperc_outbreaks_full$twomonthpd != "alltime",]

network_tplus1outbreaks <- rep(NA,num_2_ntwks)
tplus1outbreaks <- rep(NA, num_2_ntwks)

for(i in 1:num_2_ntwks){
  network_tplus1outbreaks[i] <- length(tplus1stats_full_fxd$epiunit[tplus1stats_full_fxd$inf_status == 1 & tplus1stats_full_fxd$tplus1 == (i+1)])
  tplus1outbreaks[i] <- length(outbreakepiunits[[(i+1)]])
 }

##5.27.2024: Number of Outbreaks found by all methods in t+1 months with less than x Outbreaks/more than y outbreaks

#what are reasonable x, y values?
hist(tplus1outbreaks)
quantile(tplus1outbreaks)
#    0%    25%    50%    75%   100% 
#  0.00  16.25  42.00 104.75 395.00 
#x = 15, y = 50? just to choose reasonable values

#Which t+1 months have < 15 outbreaks, >50 outbreaks?
outbreaksrare <- which(tplus1outbreaks < 15)
outbreakscommon <- which(tplus1outbreaks > 50)

#how many outbreaks did the various methods find in those months at 5% and 35% surveillance effort?
#less than 15 outbreaks, 5% found min 0% max 67% and 35% found min 33% max 100%
#more than 50 outbreaks, 5% found min 4% max 44% and 35% found min 35% max 89%


#5% - less than 15 outbreaks -> methods found max 67%
#spatial proximity
quantile(euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$twomthntwk %in% outbreaksrare & euclideandist_varr_outneighbs$pctaway == 0.05]*tplus1outbreaks[outbreaksrare], na.rm =T) #double checked that %in% preserved the order of outbreaksrare
#  0%  25%  50%  75% 100% 
#0.00 0.75 1.00 2.25 6.00 
quantile(euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$twomthntwk %in% outbreaksrare & euclideandist_varr_outneighbs$pctaway == 0.05]*100, na.rm =T)
#       0%       25%       50%       75%      100% 
# 0.000000  5.769231 14.285714 30.871212 66.666667 

#network risk factors
quantile(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$twomonthpd %in% outbreaksrare & pperc_outbreaks_only2months$P_perc == 0.05]*tplus1outbreaks[outbreaksrare], na.rm =T) 
#  0%  25%  50%  75% 100% 
#0.00 0.75 1.50 3.00 7.00 
quantile(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$twomonthpd %in% outbreaksrare & pperc_outbreaks_only2months$P_perc == 0.05]*100, na.rm =T) 
#       0%       25%       50%       75%      100% 
# 0.000000  6.818182 12.698413 25.568182 66.666667

#network proximity
quantile(wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$twomthntwk %in% outbreaksrare & wctracing_outneighbs$pctaway == 0.05]*tplus1outbreaks[outbreaksrare], na.rm =T) 
# 0%  25%  50%  75% 100% 
#0.00 1.75 3.00 4.00 5.00
quantile(wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$twomthntwk %in% outbreaksrare & wctracing_outneighbs$pctaway == 0.05]*100, na.rm =T) 
#      0%      25%      50%      75%     100% 
# 0.00000 16.66667 32.46753 46.59091 55.55556 

#35% - less than 15 outbreaks -> methods found min 33% max 100%
#spatial proximity
quantile(euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$twomthntwk %in% outbreaksrare & euclideandist_varr_outneighbs$pctaway == 0.35]*tplus1outbreaks[outbreaksrare], na.rm =T) 
#  0%  25%  50%  75% 100% 
# 1.0  4.5  6.0  9.0 11.0 
quantile(euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$twomthntwk %in% outbreaksrare & euclideandist_varr_outneighbs$pctaway == 0.35]*100, na.rm =T) 
#       0%       25%       50%       75%      100% 
# 33.33333  48.86364  59.41558  70.67308 100.00000 

#network risk factors
quantile(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$twomonthpd %in% outbreaksrare & pperc_outbreaks_only2months$P_perc == 0.35]*tplus1outbreaks[outbreaksrare], na.rm =T) 
#0%  25%  50%  75% 100% 
#1.0  5.0  7.5 10.0 12.0 
quantile(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$twomonthpd %in% outbreaksrare & pperc_outbreaks_only2months$P_perc == 0.35]*100, na.rm =T)
#0%      25%      50%      75%     100% 
#45.45455 67.46032 76.38889 83.92857 92.30769

#network proximity
quantile(wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$twomthntwk %in% outbreaksrare & wctracing_outneighbs$pctaway == 0.35]*tplus1outbreaks[outbreaksrare], na.rm =T) 
#0%  25%  50%  75% 100% 
# 1    5    6   10   12 
quantile(wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$twomthntwk %in% outbreaksrare & wctracing_outneighbs$pctaway == 0.35]*100, na.rm =T) 
# 0%      25%      50%      75%     100% 
#33.33333 54.54545 71.42857 81.81818 85.71429 

#5% - more than 50 outbreaks -> methods found min 4%, max 44%
#spatial proximity
quantile(euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$twomthntwk %in% outbreakscommon & euclideandist_varr_outneighbs$pctaway == 0.05]*tplus1outbreaks[outbreakscommon], na.rm =T) 
# 0%  25%  50%  75% 100% 
#  5   12   19   28   56 
quantile(euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$twomthntwk %in% outbreakscommon & euclideandist_varr_outneighbs$pctaway == 0.05]*100, na.rm =T) 
#0%       25%       50%       75%      100% 
#4.237288 11.580230 16.858238 21.533734 43.859649 

#network risk factors
quantile(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$twomonthpd %in% outbreakscommon & pperc_outbreaks_only2months$P_perc == 0.05]*tplus1outbreaks[outbreakscommon], na.rm =T) 
#  0%  25%  50%  75% 100% 
#   5   15   22   32   84 
quantile(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$twomonthpd %in% outbreakscommon & pperc_outbreaks_only2months$P_perc == 0.05]*100, na.rm =T) 
# 0%       25%       50%       75%      100% 
#7.894737 14.888429 19.186047 21.140303 36.923077

#network proximity
quantile(wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$twomthntwk %in% outbreakscommon & wctracing_outneighbs$pctaway == 0.05]*tplus1outbreaks[outbreakscommon], na.rm =T) 
#0%  25%  50%  75% 100% 
# 5   17   23   32   86
quantile(wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$twomthntwk %in% outbreakscommon & wctracing_outneighbs$pctaway == 0.05]*100, na.rm =T)
#0%       25%       50%       75%      100% 
#7.894737 16.666667 20.000000 22.952586 35.384615 

#35% - more than 50 outbreaks -> methods found min 35%, max 89%
#spatial proximity
quantile(euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$twomthntwk %in% outbreakscommon & euclideandist_varr_outneighbs$pctaway == 0.35]*tplus1outbreaks[outbreakscommon], na.rm =T) 
#   0%   25%   50%   75%  100% 
# 35.0  49.0  63.0 103.5 201.0 
quantile(euclideandist_varr_outneighbs$pcttplus1outbks[euclideandist_varr_outneighbs$twomthntwk %in% outbreakscommon & euclideandist_varr_outneighbs$pctaway == 0.35]*100, na.rm =T) 
#0%      25%      50%      75%     100% 
#35.64356 54.99007 61.62791 69.05983 75.38462 

#network risk factors
quantile(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$twomonthpd %in% outbreakscommon & pperc_outbreaks_only2months$P_perc == 0.35]*tplus1outbreaks[outbreakscommon], na.rm =T) 
#   0%   25%   50%   75%  100% 
# 31.0  56.0  82.0 127.5 290.0 
quantile(pperc_outbreaks_only2months$Actual_peroutbreaks[pperc_outbreaks_only2months$twomonthpd %in% outbreakscommon & pperc_outbreaks_only2months$P_perc == 0.35]*100, na.rm =T) 
#      0%      25%      50%      75%     100% 
#49.20635 64.98626 71.04072 74.62121 87.23404

#network proximity
quantile(wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$twomthntwk %in% outbreakscommon & wctracing_outneighbs$pctaway == 0.35]*tplus1outbreaks[outbreakscommon], na.rm =T) 
#0%   25%   50%   75%  100% 
#30.0  55.5  82.0 108.5 280.0 

quantile(wctracing_outneighbs$pcttplus1outbks[wctracing_outneighbs$twomthntwk %in% outbreakscommon & wctracing_outneighbs$pctaway == 0.35]*100, na.rm =T) 
#      0%      25%      50%      75%     100% 
#47.61905 61.38405 68.42105 73.40477 89.58333 