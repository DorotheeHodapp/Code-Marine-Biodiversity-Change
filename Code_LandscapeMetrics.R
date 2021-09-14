# Code to compute landscape metrics based on species distribution information provided by Aquamaps runs 
# for the different RCP scenarios

# Load packages
library(data.table)
library(raster)
library(geosphere)
library(landscapemetrics)


# Data files contain information on cell ID, species ID, species probability of occurrence, geographic information (as provided from Aquamaps runs)
df.cur # current species distribution info
df.2100.26 # future species distribution info under assumptions of RCP 2.6
df.2100.45 # future species distribution info under assumptions of RCP 4.5
df.2100.85 # future species distribution info under assumptions of RCP 8.5

#### parallel loop ####
library(doParallel)

ncores  = detectCores(logical = F)
registerDoParallel(ncores-2)

fileResults = "....csv"
fileWarning = "log_warning.txt"
fileError = "log_error.txt"


# Compute landscape metrics for current and future species distributions and RCP scenarios

set.seed(100)

start.time <- Sys.time()

result = 
  foreach(i = 1:length(speciesID), .combine=rbind, .multicombine=T, .verbose = T, .export= c('df.cur','speciesID'),                  # length(speciesID)
          .packages=c('raster', 'geosphere', 'vegan', 'data.table', 'igraph')) %dopar% {
            
            df_sp <- df.cur[df.cur$SpeciesID==speciesID[i], c("long", "lat", "Probability")]
            df_sp2100.26 <- df.2100.26[df.2100.26$SpeciesID==speciesID[i], c("long", "lat", "Probability")]
            df_sp2100.45 <- df.2100.45[df.2100.45$SpeciesID==speciesID[i], c("long", "lat", "Probability")]
            df_sp2100.85 <- df.2100.85[df.2100.85$SpeciesID==speciesID[i], c("long", "lat", "Probability")]
            
            
            tryCatch({
              
              m.cur <- suppressWarnings(rasterFromXYZ(df_sp, crs="+proj=longlat +datum=WGS84"))
              m.2100.26 <-  suppressWarnings(rasterFromXYZ(df_sp2100.26, crs="+proj=longlat +datum=WGS84"))
              m.2100.45 <-  suppressWarnings(rasterFromXYZ(df_sp2100.45, crs="+proj=longlat +datum=WGS84"))
              m.2100.85 <-  suppressWarnings(rasterFromXYZ(df_sp2100.85, crs="+proj=longlat +datum=WGS84"))
              
              
              # Compute the size of the current and future distribution area for each species (i)
              # only add up cells with prob > 0.5
              m.clean <- m.cur
              m.clean[m.clean <= 0] <-NA # remove cells with prob = 0
              
              m.clean.26 <- m.2100.26
              m.clean.26[m.clean.26 <= 0] <-NA 
              
              m.clean.45 <- m.2100.45
              m.clean.45[m.clean.45 <= 0] <-NA
              
              m.clean.85 <- m.2100.85
              m.clean.85[m.clean.85 <= 0] <-NA
              
              ### Landscape metrics
              # Current
              p_area_cv <- lsm_l_area_cv(m.clean, directions = 8)$value
              p_area_mn <- lsm_l_area_mn(m.clean, directions = 8)$value
              p_area_sd <- lsm_l_area_sd(m.clean, directions = 8)$value

              p_enn_cv <- lsm_l_enn_cv(m.clean, directions = 8)$value
              p_enn_mn <- lsm_l_enn_mn(m.clean, directions = 8)$value
              p_enn_sd <- lsm_l_enn_sd(m.clean, directions = 8)$value

              p_n <- lsm_l_np(m.clean, directions = 8)$value

              
              # 2100 RCP2.6
              p_area_cv_26 <- lsm_l_area_cv(m.clean.26, directions = 8)$value
              p_area_mn_26 <- lsm_l_area_mn(m.clean.26, directions = 8)$value
              p_area_sd_26 <- lsm_l_area_sd(m.clean.26, directions = 8)$value

              p_enn_cv_26 <- lsm_l_enn_cv(m.clean.26, directions = 8)$value
              p_enn_mn_26 <- lsm_l_enn_mn(m.clean.26, directions = 8)$value
              p_enn_sd_26 <- lsm_l_enn_sd(m.clean.26, directions = 8)$value

              p_n_26 <- lsm_l_np(m.clean.26, directions = 8)$value

              
              # 2100 RCP4.5
              p_area_cv_45 <- lsm_l_area_cv(m.clean.45, directions = 8)$value
              p_area_mn_45 <- lsm_l_area_mn(m.clean.45, directions = 8)$value
              p_area_sd_45 <- lsm_l_area_sd(m.clean.45, directions = 8)$value

              p_enn_cv_45 <- lsm_l_enn_cv(m.clean.45, directions = 8)$value
              p_enn_mn_45 <- lsm_l_enn_mn(m.clean.45, directions = 8)$value
              p_enn_sd_45 <- lsm_l_enn_sd(m.clean.45, directions = 8)$value

              p_n_45 <- lsm_l_np(m.clean.45, directions = 8)$value

              
              # 2100 RCP8.5
              p_area_cv_85 <- lsm_l_area_cv(m.clean.85, directions = 8)$value
              p_area_mn_85 <- lsm_l_area_mn(m.clean.85, directions = 8)$value
              p_area_sd_85 <- lsm_l_area_sd(m.clean.85, directions = 8)$value

              p_enn_cv_85 <- lsm_l_enn_cv(m.clean.85, directions = 8)$value
              p_enn_mn_85 <- lsm_l_enn_mn(m.clean.85, directions = 8)$value
              p_enn_sd_85 <- lsm_l_enn_sd(m.clean.85, directions = 8)$value

              p_n_85 <- lsm_l_np(m.clean.85, directions = 8)$value

              
              
              df <- data.frame(SpeciesID = speciesID[i],
                               p_area_cv = p_area_cv, p_area_mn = p_area_mn, p_area_sd = p_area_sd,
                               p_enn_cv = p_enn_cv, p_enn_mn = p_enn_mn, p_enn_sd = p_enn_sd, p_n = p_n,
                               p_area_cv_26 = p_area_cv_26, p_area_mn_26 = p_area_mn_26, p_area_sd_26 = p_area_sd_26,
                               p_enn_cv_26 = p_enn_cv_26, p_enn_mn_26 = p_enn_mn_26, p_enn_sd_26 = p_enn_sd_26, p_n_26 = p_n_26,
                               p_area_cv_45 = p_area_cv_45, p_area_mn_45 = p_area_mn_45, p_area_sd_45 = p_area_sd_45,
                               p_enn_cv_45 = p_enn_cv_45, p_enn_mn_45 = p_enn_mn_45, p_enn_sd_45 = p_enn_sd_45, p_n_45 = p_n_45,
                               p_area_cv_85 = p_area_cv_85, p_area_mn_85 = p_area_mn_85, p_area_sd_85 = p_area_sd_85,
                               p_enn_cv_85 = p_enn_cv_85, p_enn_mn_85 = p_enn_mn_85, p_enn_sd_85 = p_enn_sd_85, p_n_85 = p_n_85,)
              write.table(df, file= fileResults,
                          append=(i!=1),
                          col.names =(i==1),
                          row.names = F, sep=";")
              return(df)
              
            },error=function(e){
              msg1 = paste("ERROR on i =", paste(i), "- problem with raster for", speciesID[i])
              write.table(msg1, file= fileError, append=T, row.names = F, col.names = F)
              df = df <- data.frame(SpeciesID = speciesID[i],
                                    p_area_cv = NA, p_area_mn = NA, p_area_sd = NA,
                                    p_enn_cv = NA, p_enn_mn = NA, p_enn_sd = NA, p_n = NA,
                                    p_area_cv_26 = NA, p_area_mn_26 = NA, p_area_sd_26 = NA,
                                    p_enn_cv_26 = NA, p_enn_mn_26 = NA, p_enn_sd_26 = NA, p_n_26 = NA,
                                    p_area_cv_45 = NA, p_area_mn_45 = NA, p_area_sd_45 = NA,
                                    p_enn_cv_45 = NA, p_enn_mn_45 = NA, p_enn_sd_45 = NA, p_n_45 = NA,
                                    p_area_cv_85 = NA, p_area_mn_85 = NA, p_area_sd_85 = NA,
                                    p_enn_cv_85 = NA, p_enn_mn_85 = NA, p_enn_sd_85 = NA, p_n_85 = NA)
              write.table(df, file= fileResults, append=(i!=1),
                          col.names =(i==1),
                          row.names = F, sep=";")
              return(df)
            })
          }

