# Code to compute equatorial gaps distance for current and future species distributions

# Load packages

library(data.table)
library(raster)
library(SDMTools)
library(geosphere)
library(landscapemetrics)


#### Load data ####

# Data files contain information on cell ID, species ID, species probability of occurrence, geographic information (as provided from Aquamaps runs)
df.cur # current species distribution info
df.2100.26 # future species distribution info under assumptions of RCP 2.6
df.2100.45 # future species distribution info under assumptions of RCP 4.5
df.2100.85 # future species distribution info under assumptions of RCP 8.5

#### Compute distances between current and future distributions for each species ----

set.seed(100)

start.time <- Sys.time()

for (i in 1:length(speciesID)) {
  
  df_sp <- df.cur[df.cur$SpeciesID==speciesID[i], c("long", "lat", "Probability")]
  
  df_sp2100.26 <- df.2100.26[df.2100.26$SpeciesID==speciesID[i], c("long", "lat", "Probability")]
  df_sp2100.45 <- df.2100.45[df.2100.45$SpeciesID==speciesID[i], c("long", "lat", "Probability")]
  df_sp2100.85 <- df.2100.85[df.2100.85$SpeciesID==speciesID[i], c("long", "lat", "Probability")]
  
  
  tNlat <- quantile(df_sp$lat[df_sp$Probability >0.5 & df_sp$lat >=0], 0.01, na.rm = T) 
  tSlat <- quantile(df_sp$lat[df_sp$Probability >0.5 & df_sp$lat <=0], 0.99, na.rm = T)
  
  tNlat2050.85 <- quantile(df_sp2050.85$lat[df_sp2050.85$Probability >0.5 & df_sp2050.85$lat >=0], 0.01, na.rm = T) 
  tSlat2050.85 <- quantile(df_sp2050.85$lat[df_sp2050.85$Probability >0.5 & df_sp2050.85$lat <=0], 0.99, na.rm = T)
  
  
  tNlat2100.26 <- quantile(df_sp2100.26$lat[df_sp2100.26$Probability >0.5 & df_sp2100.26$lat >=0], 0.01, na.rm = T) 
  tSlat2100.26 <- quantile(df_sp2100.26$lat[df_sp2100.26$Probability >0.5 & df_sp2100.26$lat <=0], 0.99, na.rm = T)
  
  tNlat2100.45 <- quantile(df_sp2100.45$lat[df_sp2100.45$Probability >0.5 & df_sp2100.45$lat >=0], 0.01, na.rm = T) 
  tSlat2100.45 <- quantile(df_sp2100.45$lat[df_sp2100.45$Probability >0.5 & df_sp2100.45$lat <=0], 0.99, na.rm = T)
  
  tNlat2100.85 <- quantile(df_sp2100.85$lat[df_sp2100.85$Probability >0.5 & df_sp2100.85$lat >=0], 0.01, na.rm = T) 
  tSlat2100.85 <- quantile(df_sp2100.85$lat[df_sp2100.85$Probability >0.5 & df_sp2100.85$lat <=0], 0.99, na.rm = T)
  
  
  d_gap.cur <- distm(c(0, tNlat), c(0, tSlat), fun = distVincentyEllipsoid)/1000
  d_gap.2100.26 <- distm(c(0, tNlat2100.26), c(0, tSlat2100.26), fun = distVincentyEllipsoid)/1000
  d_gap.2100.45 <- distm(c(0, tNlat2100.45), c(0, tSlat2100.45), fun = distVincentyEllipsoid)/1000
  d_gap.2100.85 <- distm(c(0, tNlat2100.85), c(0, tSlat2100.85), fun = distVincentyEllipsoid)/1000
  d_gap.2050.85 <- distm(c(0, tNlat2050.85), c(0, tSlat2050.85), fun = distVincentyEllipsoid)/1000
  
  
  df <- data.frame(SpeciesID = speciesID[i],
                   d_gap.cur=d_gap.cur,
                   d_gap.2100.26=d_gap.2100.26,
                   d_gap.2100.45=d_gap.2100.45,
                   d_gap.2100.85=d_gap.2100.85)
  
  write.table(df, file= "....csv",
              append=(i!=1),
              col.names =(i==1),
              row.names = F, sep=";")
}




