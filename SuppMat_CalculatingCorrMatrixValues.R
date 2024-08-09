
####NETWORK CONNECTIVITY - OTHER METRICS

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

####CALCULATE OTHER NETWORK CONNECTIVITY METRICS####
#divide by the number of outbreaks observed
p_perc = c(0.05,0.1,0.15,0.20,0.25,0.3,0.35)
metrix <- c("InDegree", "InCoreness", "InStrength", "OutDegree", "OutCoreness", "OutStrength", "EigenCent", "Betweenness")
numtotalepiunits <- 54096

pperc_outbreaks_full <- data.frame(metric = rep(rep(metrix, each = length(p_perc)),(num_2_ntwks+1)),P_perc = rep(rep(p_perc,length(metrix)), (num_2_ntwks+1)), twomonthpd = rep(c(seq(1,num_2_ntwks,1), "alltime"), each = (length(p_perc)*length(metrix))), Actual_peroutbreaks = NA)


for(j in 1:num_2_ntwks){
  topOutDegreeepiunits <- topInDegreeepiunits <- topInCorenessepiunits <- topOutCorenessepiunits <- topBetweennessepiunits <- topOutStrengthepiunits <- topInStrengthepiunits <- topEigenCentepiunits <- list()
  networkmetrics_2monthsep_indiv <- networkmetrics_2monthsep[[j]] 
  #num_epiunits <- gorder(networks_2monthsep[[j]]) #removed 9.28.2023 so i could put these values on the same scale as the euclidean distance and contact tracing full network methods
  onemonthoutbreaks <- outbreakepiunits[[(j+1)]]
  
  for(i in 1:length(p_perc)){
    #which are the top P% of nodes for each metric?
    topOutDegreeepiunits[[i]] <- networkmetrics_2monthsep_indiv$Node[order(networkmetrics_2monthsep_indiv$OutDegree, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)] #topOutDegreeepiunits[[i]] <- networkmetrics_2monthsep_indiv$Node[order(networkmetrics_2monthsep_indiv$OutDegree, decreasing = TRUE)][1:(p_perc[i]*num_epiunits)]
    topInDegreeepiunits[[i]] <- networkmetrics_2monthsep_indiv$Node[order(networkmetrics_2monthsep_indiv$InDegree, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)]
    topInCorenessepiunits[[i]] <- networkmetrics_2monthsep_indiv$Node[order(networkmetrics_2monthsep_indiv$InCoreness, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)]
    topOutCorenessepiunits[[i]] <- networkmetrics_2monthsep_indiv$Node[order(networkmetrics_2monthsep_indiv$OutCoreness, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)]
    topBetweennessepiunits[[i]] <- networkmetrics_2monthsep_indiv$Node[order(networkmetrics_2monthsep_indiv$Betweenness, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)]
    topOutStrengthepiunits[[i]] <- networkmetrics_2monthsep_indiv$Node[order(networkmetrics_2monthsep_indiv$OutStrength, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)]
    topInStrengthepiunits[[i]] <- networkmetrics_2monthsep_indiv$Node[order(networkmetrics_2monthsep_indiv$InStrength, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)]
    topEigenCentepiunits[[i]] <- networkmetrics_2monthsep_indiv$Node[order(networkmetrics_2monthsep_indiv$EigenCent, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)]
    
    pperc_outbreaks_full$Actual_peroutbreaks[pperc_outbreaks_full$metric == "InDegree" & pperc_outbreaks_full$P_perc == p_perc[i] & pperc_outbreaks_full$twomonthpd == j] <- length(which(onemonthoutbreaks %in% topInDegreeepiunits[[i]]))/length(onemonthoutbreaks)
    
    pperc_outbreaks_full$Actual_peroutbreaks[pperc_outbreaks_full$metric == "InCoreness" & pperc_outbreaks_full$P_perc == p_perc[i] & pperc_outbreaks_full$twomonthpd == j] <- length(which(onemonthoutbreaks %in% topInCorenessepiunits[[i]]))/length(onemonthoutbreaks)
    
    pperc_outbreaks_full$Actual_peroutbreaks[pperc_outbreaks_full$metric == "InStrength" & pperc_outbreaks_full$P_perc == p_perc[i] & pperc_outbreaks_full$twomonthpd == j] <- length(which(onemonthoutbreaks %in% topInStrengthepiunits[[i]]))/length(onemonthoutbreaks)
    
    pperc_outbreaks_full$Actual_peroutbreaks[pperc_outbreaks_full$metric == "OutDegree" & pperc_outbreaks_full$P_perc == p_perc[i] & pperc_outbreaks_full$twomonthpd == j] <- length(which(onemonthoutbreaks %in% topOutDegreeepiunits[[i]]))/length(onemonthoutbreaks)
    
    pperc_outbreaks_full$Actual_peroutbreaks[pperc_outbreaks_full$metric == "OutCoreness" & pperc_outbreaks_full$P_perc == p_perc[i] & pperc_outbreaks_full$twomonthpd == j] <- length(which(onemonthoutbreaks %in% topOutCorenessepiunits[[i]]))/length(onemonthoutbreaks)
    
    pperc_outbreaks_full$Actual_peroutbreaks[pperc_outbreaks_full$metric == "OutStrength" & pperc_outbreaks_full$P_perc == p_perc[i] & pperc_outbreaks_full$twomonthpd == j] <- length(which(onemonthoutbreaks %in% topOutStrengthepiunits[[i]]))/length(onemonthoutbreaks)
    
    pperc_outbreaks_full$Actual_peroutbreaks[pperc_outbreaks_full$metric == "EigenCent" & pperc_outbreaks_full$P_perc == p_perc[i] & pperc_outbreaks_full$twomonthpd == j] <- length(which(onemonthoutbreaks %in% topEigenCentepiunits[[i]]))/length(onemonthoutbreaks)
    
    pperc_outbreaks_full$Actual_peroutbreaks[pperc_outbreaks_full$metric == "Betweenness" & pperc_outbreaks_full$P_perc == p_perc[i] & pperc_outbreaks_full$twomonthpd == j] <- length(which(onemonthoutbreaks %in% topBetweennessepiunits[[i]]))/length(onemonthoutbreaks)
  }
  #erase all of these, just to make sure everything is working
  topOutDegreeepiunits <- topInDegreeepiunits <- topInCorenessepiunits <- topOutCorenessepiunits <- topBetweennessepiunits <- topOutStrengthepiunits <- topInStrengthepiunits <- topEigenCentepiunits <- networkmetrics_2monthsep_indiv <- num_epiunits <- NULL
}

#want to compare that with the full network/full timeseries

#which are the top P% of nodes for each metric?
#p_perc = c(0.1,0.25,0.5,0.75)
topOutDegreeepiunits <- topInDegreeepiunits <- topInCorenessepiunits <- topOutCorenessepiunits <- topBetweennessepiunits <- topOutStrengthepiunits <- topInStrengthepiunits <- topEigenCentepiunits <- list()
for(i in 1:length(p_perc)){
  topOutDegreeepiunits[[i]] <- epiunit_networkmetrics$Node[order(epiunit_networkmetrics$OutDegree, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)]
  #changed 9.28.2023 so i could put these values on the same scale as the euclidean distance and contact tracing full network methods
  #topOutDegreeepiunits[[i]] <- epiunit_networkmetrics$Node[order(epiunit_networkmetrics$OutDegree, decreasing = TRUE)][1:(p_perc[i]*num_epiunits)] 
  topInDegreeepiunits[[i]] <- epiunit_networkmetrics$Node[order(epiunit_networkmetrics$InDegree, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)]
  topInCorenessepiunits[[i]] <- epiunit_networkmetrics$Node[order(epiunit_networkmetrics$InCoreness, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)]
  topOutCorenessepiunits[[i]] <- epiunit_networkmetrics$Node[order(epiunit_networkmetrics$OutCoreness, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)]
  topBetweennessepiunits[[i]] <- epiunit_networkmetrics$Node[order(epiunit_networkmetrics$Betweenness, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)]
  topOutStrengthepiunits[[i]] <- epiunit_networkmetrics$Node[order(epiunit_networkmetrics$OutStrength, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)]
  topInStrengthepiunits[[i]] <- epiunit_networkmetrics$Node[order(epiunit_networkmetrics$InStrength, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)]
  topEigenCentepiunits[[i]] <- epiunit_networkmetrics$Node[order(epiunit_networkmetrics$EigenCent, decreasing = TRUE)][1:(p_perc[i]*numtotalepiunits)]
}


#what percentage of the nodes in the top P% for each metric, experience outbreaks?

#make a vector containing the epiunit # of all of the epiunits that got infected
#outbreakData probably/might have duplicates in vill_ID
outbreakepiunits_nontemp <- unique(outbreakData$vill_ID)

#make that table, finally
#Actual = any epiunits that had an outbreak at ANY point
#divide by the number of outbreaks observed
metrix <- c("InDegree", "InCoreness", "InStrength", "OutDegree", "OutCoreness", "OutStrength", "EigenCent", "Betweenness")
pperc_outbreaks <- data.frame(metric = rep(metrix, each = length(p_perc)),P_perc = rep(p_perc,length(metrix)), Actual_peroutbreaks = NA)

for(i in 1:length(p_perc)){
  pperc_outbreaks$Actual_peroutbreaks[pperc_outbreaks$metric == "InDegree" & pperc_outbreaks$P_perc == p_perc[i]] <- length(which(outbreakepiunits_nontemp %in% topInDegreeepiunits[[i]]))/length(outbreakepiunits_nontemp)
  
  pperc_outbreaks$Actual_peroutbreaks[pperc_outbreaks$metric == "InCoreness" & pperc_outbreaks$P_perc == p_perc[i]] <- length(which(outbreakepiunits_nontemp %in% topInCorenessepiunits[[i]]))/length(outbreakepiunits_nontemp)
  
  pperc_outbreaks$Actual_peroutbreaks[pperc_outbreaks$metric == "InStrength" & pperc_outbreaks$P_perc == p_perc[i]] <- length(which(outbreakepiunits_nontemp %in% topInStrengthepiunits[[i]]))/length(outbreakepiunits_nontemp)
  
  pperc_outbreaks$Actual_peroutbreaks[pperc_outbreaks$metric == "OutDegree" & pperc_outbreaks$P_perc == p_perc[i]] <- length(which(outbreakepiunits_nontemp %in% topOutDegreeepiunits[[i]]))/length(outbreakepiunits_nontemp)
  
  pperc_outbreaks$Actual_peroutbreaks[pperc_outbreaks$metric == "OutCoreness" & pperc_outbreaks$P_perc == p_perc[i]] <- length(which(outbreakepiunits_nontemp %in% topOutCorenessepiunits[[i]]))/length(outbreakepiunits_nontemp)
  
  pperc_outbreaks$Actual_peroutbreaks[pperc_outbreaks$metric == "OutStrength" & pperc_outbreaks$P_perc == p_perc[i]] <- length(which(outbreakepiunits_nontemp %in% topOutStrengthepiunits[[i]]))/length(outbreakepiunits_nontemp)
  
  pperc_outbreaks$Actual_peroutbreaks[pperc_outbreaks$metric == "EigenCent" & pperc_outbreaks$P_perc == p_perc[i]] <- length(which(outbreakepiunits_nontemp %in% topEigenCentepiunits[[i]]))/length(outbreakepiunits_nontemp)
  
  pperc_outbreaks$Actual_peroutbreaks[pperc_outbreaks$metric == "Betweenness" & pperc_outbreaks$P_perc == p_perc[i]] <- length(which(outbreakepiunits_nontemp %in% topBetweennessepiunits[[i]]))/length(outbreakepiunits_nontemp)
}

#add it to the temporal networks dataframes
#pperc_outbreaks_fullwide$Actual_peroutbreaks_nontemp <- pperc_outbreaks$Actual_peroutbreaks
pperc_outbreaks_full$Actual_peroutbreaks[pperc_outbreaks_full$twomonthpd == "alltime"] <- pperc_outbreaks$Actual_peroutbreaks

#saveRDS(pperc_outbreaks_full, file = "FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthnetworks_riskfactor_fullnetworkvalues_9.28.2023.rds") #note from 10.19.2023 - the %in% statements here are fine, they just had to be reversed in the 10.13.2023 chunk because was using them to subset rows of another dataframe so it had to spit out things relevant to that dataset, the lengths (when i calculated them) were consistent regardless of the order of the %in% statements 


#Fig S2c: how does it change over time? remove the "alltime" rows first
pperc_outbreaks_only2months <- pperc_outbreaks_full[pperc_outbreaks_full$twomonthpd != "alltime",]

plot <- ggplot(data =pperc_outbreaks_only2months[pperc_outbreaks_only2months$P_perc==0.1,], aes(x=as.numeric(twomonthpd), y=Actual_peroutbreaks, group = metric))+
  geom_line(aes(color = metric))+
  ggtitle("Change in Metric Over Time")+
  xlab("2 Month Period")+
  ylab("% Outbreaks Happening in Top 10% of Nodes")+
  ylim(c(0,1)) #added 7.2.2024
#ggsave(plot, filename = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthplots/2monthnetworks_overtime_allmetrics_7.2.2024.png"), bg = "transparent", height = 10, width = 10)

#Fig S2a (2nd row): if look at the top x% of nodes, how many t+1 outbreaks do we find? avg-ing across metrics and across networks (e.g. months)
#may as well use pperc_outbreaks_full so we can say something about the full network values of this also 

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




##Fig S2d
#load in the Other Metrics 2 Month Data
#pperc_outbreaks_full_othermetrics <- readRDS(file = "FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthnetworks_riskfactor_fullnetworkvalues_9.28.2023.rds") 
pperc_outbreaks_othermetrics_only2months <- pperc_outbreaks_full_othermetrics[pperc_outbreaks_full_othermetrics$twomonthpd != "alltime",]

#load in the Degree Metric 2 Month Data
#pperc_outbreaks_full_degonly <- readRDS(file = "FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthnetworks_degreeonlyriskfactor_fullnetworkvalues_11.9.2023.rds") 
pperc_outbreaks_degonly_only2months <- pperc_outbreaks_full_degonly[pperc_outbreaks_full_degonly$twomonthpd != "alltime",]

#Put all of this into one dataframe, formatted better for calculating a correlation matrix
#I guess want to calculate one correlation matrix for each surveillance effort level
#columns = metric, rows = months (1:65)
p_perc = c(0.05,0.1,0.15,0.20,0.25,0.3,0.35)
metrix <- c("InDegree", "InCoreness", "InStrength", "OutDegree", "OutCoreness", "OutStrength", "EigenCent", "Betweenness","Degree")

corrmatrix_values <- list()
for(i in 1:length(p_perc)){
  corrmatrix_values[[i]] <- data.frame(Month = seq(1,66,1), 
                                       InDegree = pperc_outbreaks_othermetrics_only2months$Actual_peroutbreaks[pperc_outbreaks_othermetrics_only2months$metric == "InDegree" & pperc_outbreaks_othermetrics_only2months$P_perc == p_perc[i]], 
                                       InCoreness = pperc_outbreaks_othermetrics_only2months$Actual_peroutbreaks[pperc_outbreaks_othermetrics_only2months$metric == "InCoreness" & pperc_outbreaks_othermetrics_only2months$P_perc == p_perc[i]], 
                                       InStrength = pperc_outbreaks_othermetrics_only2months$Actual_peroutbreaks[pperc_outbreaks_othermetrics_only2months$metric == "InStrength" & pperc_outbreaks_othermetrics_only2months$P_perc == p_perc[i]], 
                                       OutDegree = pperc_outbreaks_othermetrics_only2months$Actual_peroutbreaks[pperc_outbreaks_othermetrics_only2months$metric == "OutDegree" & pperc_outbreaks_othermetrics_only2months$P_perc == p_perc[i]], 
                                       OutCoreness = pperc_outbreaks_othermetrics_only2months$Actual_peroutbreaks[pperc_outbreaks_othermetrics_only2months$metric == "OutCoreness" & pperc_outbreaks_othermetrics_only2months$P_perc == p_perc[i]], 
                                       OutStrength = pperc_outbreaks_othermetrics_only2months$Actual_peroutbreaks[pperc_outbreaks_othermetrics_only2months$metric == "OutStrength" & pperc_outbreaks_othermetrics_only2months$P_perc == p_perc[i]], 
                                       EigenCent = pperc_outbreaks_othermetrics_only2months$Actual_peroutbreaks[pperc_outbreaks_othermetrics_only2months$metric == "EigenCent" & pperc_outbreaks_othermetrics_only2months$P_perc == p_perc[i]], 
                                       Betweenness = pperc_outbreaks_othermetrics_only2months$Actual_peroutbreaks[pperc_outbreaks_othermetrics_only2months$metric == "Betweenness" & pperc_outbreaks_othermetrics_only2months$P_perc == p_perc[i]], 
                                       Degree = pperc_outbreaks_degonly_only2months$Actual_peroutbreaks[pperc_outbreaks_degonly_only2months$metric == "Degree" & pperc_outbreaks_degonly_only2months$P_perc == p_perc[i]])
}
#row 66 is all NaN

#Pearson or Spearman? are the values normally distributed?
hist(corrmatrix_values[[7]]$InDegree) #no
hist(corrmatrix_values[[1]]$InDegree) #no
hist(corrmatrix_values[[1]]$InStrength) #closer but no
#Spearman I guess

#Compute correlation matrices and then re-shape them for plotting
library(reshape2)
corrmatrix <- list()
melted_cormats <- list() #each of these are 81 rows long
all_corrmats <- data.frame(Var1 = NA, Var2 = NA, value = NA, P_perc = rep(c(5,10,15,20,25,30,35), each = 81))
for(i in 1:length(p_perc)){
  corrmatrix[[i]] <- cor(x = corrmatrix_values[[i]][-66,-1], method = "spearman") #removing two moonth network 66 and the month column
  melted_cormats[[i]] <- melt(corrmatrix[[i]])
  all_corrmats$Var1[all_corrmats$P_perc == p_perc[i]*100] <- melted_cormats[[i]]$Var1
  all_corrmats$Var2[all_corrmats$P_perc == p_perc[i]*100] <- melted_cormats[[i]]$Var2
  all_corrmats$value[all_corrmats$P_perc == p_perc[i]*100] <- melted_cormats[[i]]$value
}

plot <- ggplot(data = all_corrmats, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman Correlation") +
  theme_minimal()+ # minimal theme
  coord_fixed()+
  geom_text(aes(Var2, Var1, label = round(value,2)), color = "black", size = 2)+
  facet_wrap(vars(P_perc))
#ggsave(plot, filename = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthplots/RiskFactorCorrelationMatrix_overtime_allmetrics_7.2.2024.png"), bg = "transparent", height = 10, width = 10)

