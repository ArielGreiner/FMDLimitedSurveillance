##Calculating the Radius Sizes from the Spatial Proximity Data

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

####SET WORKING DIRECTORY####
setwd('/Users/akg6325/Dropbox/Github') #replace with your own directory!

#1.16.2025: Making a Plot Showing how the Radius Used for the Spatial Proximity Measurement Changes Over Time
#read in the data, we're interested in the 'radii' column from the spatial proximity method: euclideandist_varr_outneighbs
#euclideandist_varr_outneighbs <- readRDS(file = paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2months_euclidean_varr_distoutneighbourhoodmetrics_10.26.2023.rds"))

#the radii column is showing the stepsize used, in degree units
#1 degree = 111km, roughly (varies depending on how close one is to the equator)...going to just use that value because it's quite a slight variation and since each of the starting points varies in lat/long it's a bit messy
#going to only plot it until 35% even though have calculated it beyond that

p_perc = c(0.05,0.1,0.15,0.20,0.25,0.3,0.35)
pdf(paste0("FMDLimitedSurveillance/TemporalNetworks_7.2023/2monthplots/1.2025_RevisionsGraphs/SpatialProximity_RadiusPlotting_km_1.16.2025.pdf"))
plot(x = euclideandist_varr_outneighbs$twomthntwk[euclideandist_varr_outneighbs$pctaway == p_perc[1]], y = euclideandist_varr_outneighbs$radii[euclideandist_varr_outneighbs$pctaway == p_perc[1]]*111, pch = 20, col = "#d53e4f", ylim = c(0,500), xlab = "Month", ylab = "Radius (km) Searched")
lines(x = euclideandist_varr_outneighbs$twomthntwk[euclideandist_varr_outneighbs$pctaway == p_perc[1]], y = euclideandist_varr_outneighbs$radii[euclideandist_varr_outneighbs$pctaway == p_perc[1]]*111, lty = 1, col = "#d53e4f")
points(x = euclideandist_varr_outneighbs$twomthntwk[euclideandist_varr_outneighbs$pctaway == p_perc[2]], y = euclideandist_varr_outneighbs$radii[euclideandist_varr_outneighbs$pctaway == p_perc[2]]*111, pch = 20, col = "#fc8d59")
lines(x = euclideandist_varr_outneighbs$twomthntwk[euclideandist_varr_outneighbs$pctaway == p_perc[2]], y = euclideandist_varr_outneighbs$radii[euclideandist_varr_outneighbs$pctaway == p_perc[2]]*111, lty = 1, col = "#fc8d59")
points(x = euclideandist_varr_outneighbs$twomthntwk[euclideandist_varr_outneighbs$pctaway == p_perc[3]], y = euclideandist_varr_outneighbs$radii[euclideandist_varr_outneighbs$pctaway == p_perc[3]]*111, pch = 20, col = "#fee08b")
lines(x = euclideandist_varr_outneighbs$twomthntwk[euclideandist_varr_outneighbs$pctaway == p_perc[3]], y = euclideandist_varr_outneighbs$radii[euclideandist_varr_outneighbs$pctaway == p_perc[3]]*111, lty = 1, col = "#fee08b")
points(x = euclideandist_varr_outneighbs$twomthntwk[euclideandist_varr_outneighbs$pctaway == p_perc[4]], y = euclideandist_varr_outneighbs$radii[euclideandist_varr_outneighbs$pctaway == p_perc[4]]*111, pch = 20, col = "#ffffbf")
lines(x = euclideandist_varr_outneighbs$twomthntwk[euclideandist_varr_outneighbs$pctaway == p_perc[4]], y = euclideandist_varr_outneighbs$radii[euclideandist_varr_outneighbs$pctaway == p_perc[4]]*111, lty = 1, col = "#ffffbf")
points(x = euclideandist_varr_outneighbs$twomthntwk[euclideandist_varr_outneighbs$pctaway == p_perc[5]], y = euclideandist_varr_outneighbs$radii[euclideandist_varr_outneighbs$pctaway == p_perc[5]]*111, pch = 20, col = "#e6f598")
lines(x = euclideandist_varr_outneighbs$twomthntwk[euclideandist_varr_outneighbs$pctaway == p_perc[5]], y = euclideandist_varr_outneighbs$radii[euclideandist_varr_outneighbs$pctaway == p_perc[5]]*111, lty = 1, col = "#e6f598")
points(x = euclideandist_varr_outneighbs$twomthntwk[euclideandist_varr_outneighbs$pctaway == p_perc[6]], y = euclideandist_varr_outneighbs$radii[euclideandist_varr_outneighbs$pctaway == p_perc[6]]*111, pch = 20, col = "#99d594")
lines(x = euclideandist_varr_outneighbs$twomthntwk[euclideandist_varr_outneighbs$pctaway == p_perc[6]], y = euclideandist_varr_outneighbs$radii[euclideandist_varr_outneighbs$pctaway == p_perc[6]]*111, lty = 1, col = "#99d594")
points(x = euclideandist_varr_outneighbs$twomthntwk[euclideandist_varr_outneighbs$pctaway == p_perc[7]], y = euclideandist_varr_outneighbs$radii[euclideandist_varr_outneighbs$pctaway == p_perc[7]]*111, pch = 20, col = "#3288bd")
lines(x = euclideandist_varr_outneighbs$twomthntwk[euclideandist_varr_outneighbs$pctaway == p_perc[7]], y = euclideandist_varr_outneighbs$radii[euclideandist_varr_outneighbs$pctaway == p_perc[7]]*111, lty = 1, col = "#3288bd")
legend(1, 500, legend=c("5%", "10%","15%","20%", "25%","30%","35%"),col=c("#d53e4f","#fc8d59","#fee08b","#ffffbf","#e6f598","#99d594","#3288bd"), pch = rep(20,6), cex=0.8)
dev.off()

##Fig. S5b - Remaking epiunit over time plot
plot(x = seq(1,66,1), y = lengths(outbreakepiunits)[1:66], pch = 20, ylab = "# of Outbreaks", xlab = "Month")
lines(x = seq(1,66,1), y = lengths(outbreakepiunits)[1:66], lty = 1)
dev.off()