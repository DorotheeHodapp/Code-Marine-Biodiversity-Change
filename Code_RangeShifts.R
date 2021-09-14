# Code to compute distance between current and future suitable habitat


# Load packages

library(data.table)
library(raster)
library(SDMTools)
library(geosphere)
library(landscapemetrics)

# Data files contain information on cell ID, species ID, species probability of occurrence, geographic information (as provided from Aquamaps runs)
df.cur # current species distribution info
df.2100.26 # future species distribution info under assumptions of RCP 2.6
df.2100.45 # future species distribution info under assumptions of RCP 4.5
df.2100.85 # future species distribution info under assumptions of RCP 8.5


library(doParallel)

ncores  = detectCores(logical = T)
registerDoParallel(ncores)


fileResults = "....csv"
fileWarning = "log_warning.txt"
fileError = "log_error.txt"


#### Compute distances between current and future distributions for each species ----
# Distances between distribution center of gravity, leading edge and trailing edge

set.seed(100)

start.time <- Sys.time()

result = 
  foreach(i = 1:length(speciesID), .combine=rbind, .multicombine=T, .verbose = T, .export= c('df.cur','speciesID'),                  # length(speciesID)
          .packages=c('raster', 'geosphere', 'SDMTools', 'vegan', 'data.table', 'igraph')) %dopar% {
            
            df_sp <- df.cur[df.cur$SpeciesID==speciesID[i], c("long", "lat", "Probability")]
            
            df_sp2100.26 <- df.2100.26[df.2100.26$SpeciesID==speciesID[i], c("long", "lat", "Probability")]
            df_sp2100.45 <- df.2100.45[df.2100.45$SpeciesID==speciesID[i], c("long", "lat", "Probability")]
            df_sp2100.85 <- df.2100.85[df.2100.85$SpeciesID==speciesID[i], c("long", "lat", "Probability")]
            
            
            tryCatch({
              
              df_r <- suppressWarnings(rasterFromXYZ(df_sp, crs="+proj=longlat +datum=WGS84"))
             
              df_r2100.26 <-  suppressWarnings(rasterFromXYZ(df_sp2100.26, crs="+proj=longlat +datum=WGS84"))
              df_r2100.45 <-  suppressWarnings(rasterFromXYZ(df_sp2100.45, crs="+proj=longlat +datum=WGS84"))
              df_r2100.85 <-  suppressWarnings(rasterFromXYZ(df_sp2100.85, crs="+proj=longlat +datum=WGS84"))
            
              
              # Leading/trailing edges: using quantiles
              
              lNlat <- quantile(df_sp$lat[df_sp$Probability >0.5 & df_sp$lat >0], 0.95, na.rm = T) 
              tNlat <- quantile(df_sp$lat[df_sp$Probability >0.5 & df_sp$lat >0], 0.05, na.rm = T) 
              lSlat <- quantile(df_sp$lat[df_sp$Probability >0.5 & df_sp$lat <0], 0.05, na.rm = T) 
              tSlat <- quantile(df_sp$lat[df_sp$Probability >0.5 & df_sp$lat <0], 0.95, na.rm = T)
             
              lNlat2100.26 <- quantile(df_sp2100.26$lat[df_sp2100.26$Probability >0.5 & df_sp2100.26$lat >0], 0.95, na.rm = T) 
              tNlat2100.26 <- quantile(df_sp2100.26$lat[df_sp2100.26$Probability >0.5 & df_sp2100.26$lat >0], 0.05, na.rm = T) 
              lSlat2100.26 <- quantile(df_sp2100.26$lat[df_sp2100.26$Probability >0.5 & df_sp2100.26$lat <0], 0.05, na.rm = T) 
              tSlat2100.26 <- quantile(df_sp2100.26$lat[df_sp2100.26$Probability >0.5 & df_sp2100.26$lat <0], 0.95, na.rm = T)
              
              lNlat2100.45 <- quantile(df_sp2100.45$lat[df_sp2100.45$Probability >0.5 & df_sp2100.45$lat >0], 0.95, na.rm = T) 
              tNlat2100.45 <- quantile(df_sp2100.45$lat[df_sp2100.45$Probability >0.5 & df_sp2100.45$lat >0], 0.05, na.rm = T) 
              lSlat2100.45 <- quantile(df_sp2100.45$lat[df_sp2100.45$Probability >0.5 & df_sp2100.45$lat <0], 0.05, na.rm = T) 
              tSlat2100.45 <- quantile(df_sp2100.45$lat[df_sp2100.45$Probability >0.5 & df_sp2100.45$lat <0], 0.95, na.rm = T)
              
              lNlat2100.85 <- quantile(df_sp2100.85$lat[df_sp2100.85$Probability >0.5 & df_sp2100.85$lat >0], 0.95, na.rm = T) 
              tNlat2100.85 <- quantile(df_sp2100.85$lat[df_sp2100.85$Probability >0.5 & df_sp2100.85$lat >0], 0.05, na.rm = T) 
              lSlat2100.85 <- quantile(df_sp2100.85$lat[df_sp2100.85$Probability >0.5 & df_sp2100.85$lat <0], 0.05, na.rm = T) 
              tSlat2100.85 <- quantile(df_sp2100.85$lat[df_sp2100.85$Probability >0.5 & df_sp2100.85$lat <0], 0.95, na.rm = T)

              
              d_lN.2100.26 <- distm(c(0, lNlat), c(0, lNlat2100.26), fun = distVincentyEllipsoid)/1000
              d_lN.2100.26 <- ifelse(lNlat2100.26 >= lNlat, d_lN.2100.26 , -d_lN.2100.26) 
              d_tN.2100.26 <- distm(c(0, tNlat), c(0, tNlat2100.26), fun = distVincentyEllipsoid)/1000
              d_tN.2100.26 <- ifelse(tNlat2100.26 > tNlat, -d_tN.2100.26, d_tN.2100.26)
              d_lS.2100.26 <- distm(c(0, lSlat), c(0, lSlat2100.26), fun = distVincentyEllipsoid)/1000
              d_lS.2100.26 <- ifelse(lSlat2100.26 <= lSlat, d_lS.2100.26, -d_lS.2100.26)
              d_tS.2100.26 <- distm(c(0, tSlat), c(0, tSlat2100.26), fun = distVincentyEllipsoid)/1000
              d_tS.2100.26 <- ifelse(tSlat2100.26 < tSlat, -d_tS.2100.26, d_tS.2100.26)
              
              d_lN.2100.45 <- distm(c(0, lNlat), c(0, lNlat2100.45), fun = distVincentyEllipsoid)/1000
              d_lN.2100.45 <- ifelse(lNlat2100.45 >= lNlat, d_lN.2100.45 , -d_lN.2100.45) 
              d_tN.2100.45 <- distm(c(0, tNlat), c(0, tNlat2100.45), fun = distVincentyEllipsoid)/1000
              d_tN.2100.45 <- ifelse(tNlat2100.45 > tNlat, -d_tN.2100.45, d_tN.2100.45)
              d_lS.2100.45 <- distm(c(0, lSlat), c(0, lSlat2100.45), fun = distVincentyEllipsoid)/1000
              d_lS.2100.45 <- ifelse(lSlat2100.45 <= lSlat, d_lS.2100.45, -d_lS.2100.45)
              d_tS.2100.45 <- distm(c(0, tSlat), c(0, tSlat2100.45), fun = distVincentyEllipsoid)/1000
              d_tS.2100.45 <- ifelse(tSlat2100.45 < tSlat, -d_tS.2100.45, d_tS.2100.45)
              
              d_lN.2100.85 <- distm(c(0, lNlat), c(0, lNlat2100.85), fun = distVincentyEllipsoid)/1000
              d_lN.2100.85 <- ifelse(lNlat2100.85 >= lNlat, d_lN.2100.85 , -d_lN.2100.85) 
              d_tN.2100.85 <- distm(c(0, tNlat), c(0, tNlat2100.85), fun = distVincentyEllipsoid)/1000
              d_tN.2100.85 <- ifelse(tNlat2100.85 > tNlat, -d_tN.2100.85, d_tN.2100.85)
              d_lS.2100.85 <- distm(c(0, lSlat), c(0, lSlat2100.85), fun = distVincentyEllipsoid)/1000
              d_lS.2100.85 <- ifelse(lSlat2100.85 <= lSlat, d_lS.2100.85, -d_lS.2100.85)
              d_tS.2100.85 <- distm(c(0, tSlat), c(0, tSlat2100.85), fun = distVincentyEllipsoid)/1000
              d_tS.2100.85 <- ifelse(tSlat2100.85 < tSlat, -d_tS.2100.85, d_tS.2100.85)
            
              ### Suitable habitat amount
              
              m.cur <- round(df_r) # keep cells > 0.5 (scores range from 0 to 1)
              m.2100.26 <- round(df_r2100.26)
              m.2100.45 <- round(df_r2100.45)
              m.2100.85 <- round(df_r2100.85)
              
              
              # Compute the size of the current and future distribution area for each species (i)
              # only add up cells with prob > 0.5
              
              m.clean <- m.cur
              m.clean[m.clean <= 0] <-NA # remove cells with prob = 0
              cell.size.cur <- area(m.clean, na.rm = T)
              cell.size.cur <- cell.size.cur[!is.na(cell.size.cur)]
              a.size.cur <- sum(cell.size.cur, na.rm = T) # in km2
              
              m.clean.2050.85 <- m.2050.85
              m.clean.2050.85[m.clean.2050.85 <= 0] <-NA
              cell.size.2050.85 <- area(m.clean.2050.85, na.rm = T)
              cell.size.2050.85 <- cell.size.2050.85[!is.na(cell.size.2050.85)]
              a.size.2050.85 <- sum(cell.size.2050.85, na.rm = T) # in km2
              
              m.clean.26 <- m.2100.26
              m.clean.26[m.clean.26 <= 0] <-NA 
              cell.size.2100.26 <- area(m.clean.26, na.rm = T)
              cell.size.2100.26 <- cell.size.2100.26[!is.na(cell.size.2100.26)]
              a.size.2100.26 <- sum(cell.size.2100.26, na.rm = T) # in km2
              
              m.clean.45 <- m.2100.45
              m.clean.45[m.clean.45 <= 0] <-NA
              cell.size.2100.45 <- area(m.clean.45, na.rm = T)
              cell.size.2100.45 <- cell.size.2100.45[!is.na(cell.size.2100.45)]
              a.size.2100.45 <- sum(cell.size.2100.45, na.rm = T) # in km2
              
              m.clean.85 <- m.2100.85
              m.clean.85[m.clean.85 <= 0] <-NA
              cell.size.2100.85 <- area(m.clean.85, na.rm = T)
              cell.size.2100.85 <- cell.size.2100.85[!is.na(cell.size.2100.85)]
              a.size.2100.85 <- sum(cell.size.2100.85, na.rm = T) # in km2

              
                            df <- data.frame(SpeciesID = speciesID[i],
                             d_lN.2100.26 = d_lN.2100.26, d_lN.2100.45 = d_lN.2100.45, d_lN.2100.85 = d_lN.2100.85,
                               d_tN.2100.26 = d_tN.2100.26, d_tN.2100.45 = d_tN.2100.45, d_tN.2100.85 = d_tN.2100.85,
                               d_lS.2100.26 = d_lS.2100.26, d_lS.2100.45 = d_lS.2100.45, d_lS.2100.85 = d_lS.2100.85,
                               d_tS.2100.26 = d_tS.2100.26, d_tS.2100.45 = d_tS.2100.45, d_tS.2100.85 = d_tS.2100.85,
                               sp.area.cur = a.size.cur, sp.area.2100.26 = a.size.2100.26, sp.area.2100.45 = a.size.2100.45, sp.area.2100.85=a.size.2100.85)
              write.table(df, file= fileResults,
                          append=(i!=1),
                          col.names =(i==1),
                          row.names = F, sep=";")
              return(df)
              
            },error=function(e){
              msg1 = paste("ERROR on i =", paste(i), "- problem with raster for", speciesID[i])
              write.table(msg1, file= fileError, append=T, row.names = F, col.names = F)
              df = df <- data.frame(SpeciesID = speciesID[i],
                                    d_lN.2100.26 = NA, d_lN.2100.45 = NA, d_lN.2100.85 = NA,
                                    d_tN.2100.26 = NA, d_tN.2100.45 = NA, d_tN.2100.85 = NA,
                                    d_lS.2100.26 = NA, d_lS.2100.45 = NA, d_lS.2100.85 = NA,
                                    d_tS.2100.26 = NA, d_tS.2100.45 = NA, d_tS.2100.85 = NA,
                                    sp.area.cur = NA, sp.area.2100.26 = NA, sp.area.2100.45 = NA, sp.area.2100.85=NA)
              write.table(df, file= fileResults, append=(i!=1),
                          col.names =(i==1),
                          row.names = F, sep=";")
              return(df)
            })
          }

