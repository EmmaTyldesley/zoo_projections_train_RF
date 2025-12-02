# check stability of error and predictor importance:
# train rf 10 times with different starting seeds
# (time consuming: use data subset initially for speed, then run over weekend)

# -----------------------------------------
# EJT December 2025 - tidied up for manuscript submission

# required packages
library(randomForest)
library(dplyr) 
library(lubridate)

rm(list = ls())

# make results repeatable - reset seed here in case one model affects next
#set.seed(9)
n_reps <- 10

# stripped down workflow

# ---- load cleaned ZED and AMM7 dataset ---
load("data/logzed_and_predictors.RData")
dataset <- logzed_with_amm7
n_samples <- nrow(dataset)

# ---- define zooplankton classes ----
zooClasses <- c("small","Cfin","Chel")
n_zooClasses <- length(zooClasses)
opt_mtry_Z <- c(3,7,8)

# define models
mdl_name <- c("NEAtlantic") 
# subsample for full model
doSubset <- FALSE
subset_size <- 0.5
if (doSubset & mdl_name=="NEAtlantic") {
  subsample_rows <- sample(1:n_samples,subset_size*n_samples) # get random sub-sample of total obs
  datasubset <- dataset[subsample_rows,] 
} else {
  datasubset <- dataset
}

rf_string <- "full_with_no3"
datasubset <- datasubset[datasubset$sss>30,] # take out low salinities - is this too low? could do with a buffer
datasubset <- datasubset[datasubset$no3_winter<20,] #take out high nitrate

print(paste("datasubset has",nrow(datasubset),"records"))

for (j in 3) { #1:n_zooClasses) {
  
  print(paste("zoo class:",zooClasses[j]))
  datasubset$Z <- datasubset[[zooClasses[j]]] # this could be log or not
  
  # train model n_reps times
  
  OOB_error <- rep(NA, n_reps)
  var_explained <- rep(NA,n_reps)
  rsq <- rep(NA,n_reps)
  
  predictors <- sort(c("yday","sst","mo_sst_anom","sss","logchl","logdepth","SPG_proxy_anom_annual","logmld","bloom_start","bloom_duration","no3_winter","no3_mo_anom"))
  n_p <- length(predictors)
  blank<-rep(NA,n_p)
  imp_all <-data.frame(predictors=predictors,
                       imp_rep1=blank,
                       imp_rep2=blank,
                       imp_rep3=blank,
                       imp_rep4=blank,
                       imp_rep4=blank)
  
  for (i in 1:n_reps) {
    # fit model
    rf <- randomForest(Z ~ yday + sst + mo_sst_anom + sss + logchl + logdepth + SPG_proxy_anom_annual + logmld + bloom_start + bloom_duration + no3_winter + no3_mo_anom,
                       data=datasubset,
                       importance=TRUE,
                       ntree=500,
                       corr.bias=FALSE,
                       mtry=opt_mtry_Z[j])
    print(rf)
    
    # get predictions and observations
    pred <- predict(rf,datasubset)
    obs <- datasubset$Z
    
    # plot observations vs predictions
    #plot(obs,pred,main=zooClasses[j])
    
    #get mean square error (residuals)
    OOB_error[i] <- rf$mse[length(rf$mse)]
    var_explained[i] <- rf$rsq[length(rf$rsq)]*100
    rsq[i]<-cor(obs,pred)^2 # rsquared: better metric on the testing data
    
    # plot variable importance
    varImpPlot(rf,
               type=1,
               scale=TRUE,
               cex=0.8,
               main=paste(zooClasses[j],"rep ",i))
    
    imp <-as.data.frame(importance(rf))
    imp <- imp[order(row.names(imp)), ]
    imp_all[,(i+1)] <- imp$`%IncMSE`
    
  }
  #save results
  write.csv(imp_all,
            file=paste0("sensitivity analysis/","imp_",zooClasses[j]),
            row.names = FALSE)
  write.csv(OOB_error,
            file=paste0("sensitivity analysis/","error_",zooClasses[j]),
            row.names = FALSE)
  write.csv(var_explained,
            file=paste0("sensitivity analysis/","var_expl_",zooClasses[j]),
            row.names = FALSE)
  write.csv(rsq,
            file=paste0("sensitivity analysis/","rsq_",zooClasses[j]),
            row.names = FALSE)
}
