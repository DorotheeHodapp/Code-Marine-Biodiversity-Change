## Random Forest analyses assessing the impact of envionrmnetal variables and their changes on specie socmpositional turnover 
## between the beginning and the end of 21st century
## Code adapted from https://uc-r.github.io/random_forests

# Load packages

library(data.table)  # for fread 
library(dplyr)
library(ggplot2)
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many machine learning models
library(raster)
library(vip)          # for extracting variable importance information from random forest object
library(pdp)          # for extracting variable importance information from random forest object
library(doParallel)   # load the parallel backend

#### Load data ----
## Data file: data_rf_ls (Data frame containing cell ID, environmental variables, response variable (turnover))

#### Full grid search using package "ranger" ----

# Loop over all reponse variables 
model_OOB_ls <- list()

# hyperparameter grid for search
hyper_grid <- expand.grid(
  mtry       = seq(4, 10, by = 2),   # mtry: the number of variables to randomly sample as candidates at each split
  node_size  = seq(3, 9, by = 2),     # node size: minimum number of samples within the terminal nodes
  sample_size = c(.6, .7, .8),    # sample size: the number of samples to train on
  OOB_RMSE   = 0
)

# total number of combinations
nrow(hyper_grid)

start.time <- Sys.time()   

for(j in 1:length(data_rf_ls)){
  
  data_rf <- data_rf_ls[[j]]   ## choose response variable
  response <- response_var[j]
  
  for(i in 1:nrow(hyper_grid)) {       
    
    # train model
    model <- ranger(
      formula         = paste(response," ~ ."),      #paste(response," ~ .")
      data            = data_rf, 
      num.trees       = 1000,
      mtry            = hyper_grid$mtry[i],
      importance      = "permutation", 
      min.node.size   = hyper_grid$node_size[i],
      sample.fraction = hyper_grid$sample_size[i],
      seed            = 123
    )
    
    # add OOB error to grid
    hyper_grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
    model_OOB_ls[[j]] <- hyper_grid
  }
  
}


#### Rerun analysis for best performing parameter set ----

model_ls <- list()

start.time <- Sys.time()  
for (j in 1:length(data_rf_ls)){
  
  data_rf <- data_rf_ls[[j]]   ## choose response variable
  response <- response_var[j]
  
  ranking <- model_OOB_ls[[j]] %>% 
    dplyr::arrange(OOB_RMSE) %>%
    head(10)
  
  
  
  OOB_RMSE <- vector(mode = "numeric", length = 10)
  
  for(i in seq_along(OOB_RMSE)) {
    print(i)
    
    optimal_ranger <- ranger(
      formula         = paste(response," ~ ."), 
      data            = data_rf, 
      num.trees       = 1000,           
      mtry            = ranking$mtry[1],
      min.node.size   = ranking$node_size[1],
      sample.fraction = ranking$sample_size[1],
      importance      = 'permutation',
      write.forest    = TRUE,
      keep.inbag      = TRUE
    )
    
    OOB_RMSE[i] <- sqrt(optimal_ranger$prediction.error)
    model_ls[[j]] <- optimal_ranger
  }
  
  hist(OOB_RMSE, breaks = 20)
} 

# Variable imoprtance and model performance
model_ls[[1]]$r.squared 
sqrt(model_ls[[1]]$prediction.error)

# Turnover
var.imp <- as.data.frame(unlist(model_ls[[1]]$variable.importance)) 
names(var.imp) <- "importance"

png(paste("rf_varimp_", response_var[1],"_",year,"_",rcp,"_only_delta_env_var.png", sep = ""), width = 2000, height = 2000, res = 300)
ggplot(var.imp, aes(reorder(rownames(var.imp),importance),importance)) +
  geom_col() +
  coord_flip() + 
  ylab("Variable importance") + 
  xlab("Environmental variable") + 
  ggtitle("Species turnover (beginning to end 21th century) RCP2.6")
dev.off()

# Variable importance plots using package "pdp"

partial_varimp_dSST_dPP = partial(model_ls[[1]],
                                  pred.var = c('d_SST', 'd_PrimProd'), grid.resolution = 50,
                                  progress= "text", parallel=F)

partial_varimp_dSST_dIceCon = partial(model_ls[[1]],
                                      pred.var = c('d_SST', 'd_IceCon'), grid.resolution = 50,
                                      progress= "text", parallel=F)

partial_varimp_dPP_dSal = partial(model_ls[[1]],
                                  pred.var = c('d_PrimProd', 'd_Sal'), grid.resolution = 50,
                                  progress= "text", parallel=F)


png(paste("var_imp_PP_dPP_",rcp,".png", sep = ""), width = 2000, height = 2000, res = 300)
plotPartial(partial_varimp_PP_dPP,
            xlab = list(cex = 2),
            ylab = list(cex = 2))
dev.off()

png(paste("var_imp_SST_PP_",rcp,".png", sep = ""), width = 2000, height = 2000, res = 300)
plotPartial(partial_varimp_SST_PP, 
            xlab = list(cex= 2),
            ylab = list(cex= 2))
dev.off()

png(paste("var_imp_SST_SeaIce_",rcp,".png", sep = ""), width = 2000, height = 2000, res = 300)
plotPartial(partial_varimp_SST_IceCon, 
            xlab = list(cex= 2),
            ylab = list(cex= 2))
dev.off()