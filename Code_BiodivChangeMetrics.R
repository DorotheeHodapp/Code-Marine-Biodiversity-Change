## Calculation of diversity change metrics between beginnign and end of 21st century based on Aquamaps species distribution data 
## for three different CO2 emission scenarios (RCP 2.6, RCP 4.5, RCP 8.5)

## Calculate diversity change metrics (current - future) for global marine species inventory (Aquamaps output)
## Metrics: current species richness, future species richness, future species loss, future species gain, Jaccard turnover,
##          nestedness turnover, replacement turnover (sensu Baselga et al. 2010), relative species loss 

# Load packages
library(data.table)
library(betapart)
library(stringr)
library(dplyr)
library(raster)
library(foreach)
library(doMC)
library(tidyr)
library(ggplot2)

# Determine number of cores
ncores = detectCores(logical = T)
registerDoMC(cores=ncores)



# Set year to be used in output file name
future.year  <- "2100"
# Set RCP scenario (RCP26, RCP45, RCP85)
rcp <- "..."
# set probability threshold to be used
sel.prob    <- 0.5
# Set minimum number of species per cell
min.spec    <- 5
# Set what subset of species to analyse
spec.type <- "all"


# Load data
HSPEC <- fread("....csv", header=TRUE) # current species distribution containing species ID, cell ID, coordinates, probability of occurrence
HSPECfut  <- fread("....csv",header=TRUE) # future species distribution containing species ID, cell ID, coordinates, probability of occurrence
HCAF <- fread("....csv", header=TRUE) # half-degree authority file (env param info) in 0.5° resolution


HSPEC_05 <- HSPEC[HSPEC$Probability >= sel.prob,] # Only select species with current occurrence probability > sel.prob
HSPECfut_05 <- HSPECfut[HSPECfut$Probability >= sel.prob,]# Only select species with future occurrence probability > sel.prob

#### parallel loop ####

start_time <- Sys.time()

# Process one cell at a time
result =
  foreach(i = 1 : length(HCAF$LOICZID),.combine=rbind, .multicombine=TRUE, .verbose = T) %dopar% {   
    
    CellID        <- HCAF$LOICZID[i] # get cellID
    Csquare       <- HCAF$CsquareCode[i] # get Csquarecode of this cell
    
    # Process current occurrences 
    # get vector with SpeciesIDs occuring in this cell with selected habitat suitability threshold
    cell.cur        <- HSPEC_05[HSPEC_05$CsquareCode == Csquare,]
    cell.fut        <- HSPECfut_05[HSPECfut_05$CsquareCode == Csquare,]
    
    # skip cells with less than min.spec species for current and future conditions
    if(((nrow(cell.cur) < min.spec) & (nrow(cell.fut) < min.spec)) == T) {output <- data.frame()}
    else {
      
      # get counts of species for fish and other species for current year
      spec.cur         <- nrow(cell.cur)
      
      # calculate number of species lost or newly arrived
      spec.fut       <- nrow(cell.fut)  # get number of species in future
      fut.loss       <- spec.cur - length(intersect(cell.cur$SpeciesID,cell.fut$SpeciesID)) # species loss
      fut.new        <- spec.fut - spec.cur + fut.loss  # immigration
      
      # calculate relative proportion of species lost 
      rel.spec.loss <- fut.loss/spec.cur 
      
      # use presence/absence based turnover partitioning (betapart package, Baselga et al. 2010) to distinguish between species loss and replacement
      # join current and future species information to get overall species pool for this cell
      spec.comb <- full_join(cell.cur, cell.fut, by = "SpeciesID")
      spec.comb.p.a <- transmute(spec.comb, p.a.cur = ifelse(is.na(CsquareCode.x), 0, 1), p.a.fut = ifelse(is.na(CsquareCode.y), 0, 1))
      
      # calculate turnover components (output: beta.jtu, beta.jne, beta.jac)
      to <- beta.temp(as.matrix(t(spec.comb.p.a$p.a.cur)),as.matrix(t(spec.comb.p.a$p.a.fut)), index.family = "jaccard")
      
      # append data to outfile
      output = data.frame(CellID, Csquare,
                          HCAF$CenterLat[HCAF$CsquareCode==Csquare],
                          HCAF$CenterLong[HCAF$CsquareCode==Csquare],
                          HCAF$LME[HCAF$CsquareCode==Csquare],
                          spec.cur, spec.fut, fut.loss,
                          fut.new, rel.spec.loss, to$beta.jtu, to$beta.jne, to$beta.jac)
    }
    return(output) 
  }

write.table(result, file=outfile, row.names = F)  

