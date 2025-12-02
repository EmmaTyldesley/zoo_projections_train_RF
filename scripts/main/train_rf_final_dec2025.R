# --- train RF models - final iteration ---
# -----------------------------------------
# EJT December 2025 - tidied up for manuscript submission

# model: NE Atlantic
# predictors:
#   depth
#   daily sst, sss, chl, diato frac, dino frac, mixed layer depth
#   monthly sst anomaly & NO3 anomaly
#   annual SPG proxy
#   annual bloom start and duration
#   annual winter NO3
# Note: exclude bloom phenology for HADGEM projection "reduced model" since daily P dynamics not available
# steps:
#   load zed and predictors dataset; zed is log10(zed+offset), where offset=0.5(min(non zero ZED))
#   these offsets need to be loaded when running the forward models
#   fit rf model for each zooplankton class: aggregated small copepods, Calanus finmarchicus and C. helgolandicus

# required packages
library(randomForest)
library(dplyr) 
library(lubridate)

# clear workspace
rm(list = ls())
# make results repeatable
set.seed(9)

saveData <- TRUE

# ---- load cleaned ZED and AMM7 dataset (see pre_processing/clean_up_data.R)---
load("data/logzed_and_predictors.RData")
dataset <- logzed_with_amm7
n_samples <- nrow(dataset)
print("loaded logzed dataset")

# define zooplankton classes
zooClasses <- c("small","Cfin","Chel")
n_zooClasses <- length(zooClasses)
opt_mtry_Z <- c(3,7,8) # optimal mtry as found by varying from 1:n_predictors to minimize OOB error (see sensitivity_analysis/sensitivity_mtry.R)

# define model
mdl_name <- "NEAtlantic"

# sub-sample if needed for training
doSubset <- FALSE
subset_size <- 0.1
if (doSubset) {
  subsample_rows <- sample(1:n_samples,subset_size*n_samples) # get sub-sample of total obs, 
  datasubset <- dataset[subsample_rows,] 
}

# exclude low salinities (<30)
print("Excluding coastal salinities")
print("--------")
datasubset <- datasubset[datasubset$sss>30,]

# exclude high NO3 (>20)
print("Excluding coastal nitrate")
print("--------")
datasubset <- datasubset[datasubset$no3_winter<20,]

print(paste("datasubset has",nrow(datasubset),"records"))
print("--------")

# train rf models
rf_string <- c("full_with_no3", # (for RECICLE projections) 
               "reduced_with_no3") # (for HADGEM projections)

for (r in 1:length(rf_string)) { # loop over models
  
  print(paste("Using RF model ",rf_string[r]))
  print("--------")
  
  # bias correction:
  # initialise m_ROE and b_ROE vectors
  # ROE = regression of obs on estimated values (Belitz & Stackelberg (2021). https://doi.org/10.1016/j.envsoft.2021.105006)
  m_ROE <- double(length=n_zooClasses)
  b_ROE <- double(length=n_zooClasses)
  
  # j loops over zooplankton groups
  for (j in 1:n_zooClasses) { 
    
    # which Z group
    print("--------")
    print(paste("zoo class:",zooClasses[j]))
    datasubset$Z <- datasubset[[zooClasses[j]]] # this is log10+offset
    
    # fit model
    if (rf_string[r]=="full_with_no3") { 
      rf <- randomForest(Z ~ yday + sst + mo_sst_anom + sss + logchl + logdepth + SPG_proxy_anom_annual + logmld + bloom_start + bloom_duration + no3_winter + no3_mo_anom,
                         data=datasubset,
                         importance=TRUE,
                         ntree=500,
                         mtry=opt_mtry_Z[j],
                         corr.bias=FALSE)
      
    }  else if (rf_string[r]=="reduced_with_no3") { # for HADGEM
      rf <- randomForest(Z ~ yday + sst + mo_sst_anom + sss + logchl + logdepth + SPG_proxy_anom_annual + logmld  + no3_winter + no3_mo_anom,
                         data=datasubset,
                         importance=TRUE,
                         ntree=500,
                         mtry=opt_mtry_Z[j],
                         corr.bias=FALSE)
    }
    
    print(rf)
    
    # save model
    rf_name = paste0('output/rf_',rf_string[r],'_',mdl_name,'_',zooClasses[j],'.Rdata')
    print(paste("model name:",rf_name[r]))
    if (saveData) {
      save(rf,file=rf_name)
    }
    
    # get predictions and observations
    pred <- predict(rf,datasubset)
    obs <- datasubset$Z
    
    # plot observations vs predictions
    plot(obs,pred,main=zooClasses[j])
    
    # plot variable importance
    varImpPlot(rf,
               type=1,
               scale=TRUE,
               cex=0.8,
               main='')
    
    # get importance data
    imp<-as.data.frame(importance(rf,type=1,scale=TRUE))
    predictors <-rownames(imp) # names of predictors
    n_predictors <- length(predictors)
    
    # get partial dependencies
    print('---getting partial dependencies:---')
    for (k in 1:n_predictors) {
      print(predictors[k])
      pd <- partialPlot(rf,datasubset,predictors[k],plot=FALSE)
      savename <- paste0("output/pd_",rf_string[r],'_',mdl_name,"_",zooClasses[j],"_",predictors[k],".csv")
      if (saveData) { 
      write.csv(pd, file = savename,row.names=FALSE)
        }
    }
    
    # Bias correction:
    # note: their Y^TR_ML is our pred; their Y^TR_OBS is our obs
    m_ROE[j] <- cov(pred,obs)/var(pred)
    b_ROE[j] <- mean(obs) - m_ROE[j]*mean(pred)
    
  }
  
  # save bias correction values for applying to predictions
  ROE_values <- data.frame(zooClasses,m_ROE,b_ROE)
  dataname <- paste0("data/ROE_values",mdl_name,"_",rf_string[r],".csv")
  if (saveData) { 
    write.csv(ROE_values,file=dataname,row.names = FALSE)
  }
}

