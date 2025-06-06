####CODE TO COMPARE THE NETWORK CHARACTERISTICS OF THE 2 MONTH NETWORKS AND ALL-TIME NETWORK

####SET WORKING DIRECTORY####
setwd('/Users/akg6325/Dropbox/Github') #replace with your own directory!

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

#####CHARACTERISTICS OF ALL-TIME NETWORK
#epiunit_networkmetrics #49580 nodes
median(epiunit_networkmetrics$Degree) #105
as.numeric(quantile(epiunit_networkmetrics$Degree)[4]) - as.numeric(quantile(epiunit_networkmetrics$Degree)[2]) #187
quantile(epiunit_networkmetrics$Degree) #25%: 39, 75%: 226 (226-39 = 187)
median(epiunit_networkmetrics$Betweenness) #12640.02
quantile(epiunit_networkmetrics$Betweenness) #25%: 2.012e+03, 75%: 5.08e+04 
#as.numeric(quantile(epiunit_networkmetrics$Betweenness)[4]) - as.numeric(quantile(epiunit_networkmetrics$Betweenness)[2]) #48751.82

#####AVG CHARACTERISTICS OF ALL 2 MONTH NETWORKs
#networkmetrics_2monthsep #list of 66 elements, each of which takes the same structure as epiunit_networkmetrics
avgdegree <- avgbetweenness <- rep(NA,66)
upperdegree <- upperbetweenness <- rep(NA,66)
lowerdegree <- lowerbetweenness <- rep(NA,66)
#iqrdegree <- iqrbetweenness <- rep(NA,66)
sizes <- rep(NA,66)
for(i in 1:66){
  sizes[i] <- dim(networkmetrics_2monthsep[[i]])[1]
  avgdegree[i] <- median(networkmetrics_2monthsep[[i]]$Degree)
  avgbetweenness[i] <- median(networkmetrics_2monthsep[[i]]$Betweenness)
  lowerdegree[i] <- as.numeric(quantile(networkmetrics_2monthsep[[i]]$Degree)[2])
  upperdegree[i] <- as.numeric(quantile(networkmetrics_2monthsep[[i]]$Degree)[4])
  lowerbetweenness[i] <- as.numeric(quantile(networkmetrics_2monthsep[[i]]$Betweenness)[2])
  upperbetweenness[i] <- as.numeric(quantile(networkmetrics_2monthsep[[i]]$Betweenness)[4])
  #iqrdegree[i] <- as.numeric(quantile(networkmetrics_2monthsep[[i]]$Degree)[4]) - as.numeric(quantile(networkmetrics_2monthsep[[i]]$Degree)[2])
  #iqrbetweenness[i] <- as.numeric(quantile(networkmetrics_2monthsep[[i]]$Betweenness)[4]) - as.numeric(quantile(networkmetrics_2monthsep[[i]]$Betweenness)[2])
}
mean(sizes) #37744.71
#quantile(sizes) #25%: 36212.25, 75%: 39290.25
mean(avgdegree) #6.08
mean(avgbetweenness) #2078.90
mean(lowerbetweenness) #0
mean(upperbetweenness) #39,295.8
mean(upperdegree) #14.35
mean(lowerdegree) #2.47
#mean(iqrdegree) #11.88
#mean(iqrbetweenness) #39295.8







