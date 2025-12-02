# check ability to extrapolate in time
# including stability of predictor importance ranking
# -----------------------------------------
# EJT December 2025 - tidied up for manuscript submission

# required packages
library(randomForest)
library(dplyr) 
library(lubridate)
library(Metrics) # for rmse

rm(list = ls())

# ---- load cleaned ZED and AMM7 dataset ---
load("data/logzed_and_predictors.RData")
dataset <- logzed_with_amm7
n_samples <- nrow(dataset)

# ---- define zooplankton classes ----
zooClasses <- c("small","Cfin","Chel")
n_zooClasses <- length(zooClasses)
opt_mtry_Z <- c(3,7,8) # optimal mtry


# define models
mdl_name <- c("NEAtlantic") 
# subsample for full model
doSubset <-  TRUE
subset_size <- 1
if (doSubset & mdl_name=="NEAtlantic") {
  subsample_rows <- sample(1:n_samples,subset_size*n_samples) # get random sub-sample of total obs
  datasubset <- dataset[subsample_rows,] 
} else {
  datasubset <- dataset
}

print(paste("datasubset has",nrow(datasubset),"records"))

rf_string <- "full_with_no3"
datasubset <- datasubset[datasubset$sss>28,] # take out low salinities

# --- define temporal blocks - try 7 year blocks
# note SPG variability is longer than this so won't be useful for that variable
blocks <- list(
  block1=data.frame(ymin=1998, ymax=2005),
  block2=data.frame(ymin=2006, ymax=2012),
  block3=data.frame(ymin=2013, ymax=2019) )

n_blocks <- length(blocks)
blocks_col <-rainbow(n_blocks) # define plotting colours


# --- initialize metrics table
output_metrics <- data.frame(block_num=1:n_blocks,
                             OOB_error_training=rep(NA, n_blocks),
                             rOOB_error_training=rep(NA,n_blocks),
                             var_explained_training=rep(NA,n_blocks),
                             rmse_training=rep(NA,n_blocks),
                             rsq_training=rep(NA,n_blocks),
                             rmse_test=rep(NA,n_blocks),
                             rsq_test=rep(NA,n_blocks))

output<-list(small=output_metrics,
             Cfin=output_metrics,
             Chel=output_metrics)

# initialise importance output
predictors <- sort(c("yday","sst","mo_sst_anom","sss","logchl","logdepth","SPG_proxy_anom_annual","logmld","bloom_start","bloom_duration","no3_winter","no3_mo_anom"))
n_p <- length(predictors)
blank<-rep(NA,n_p)
imp_all <-data.frame(predictors=predictors,
                     imp_block1=blank,
                     imp_block2=blank,
                     imp_block3=blank)

output_imp <-list(small=imp_all,
             Cfin=imp_all,
             Chel=imp_all)

# --- loop over blocks
for (b in 1:n_blocks) {
  
  print(paste("---- block",b,"----"))
  
  # put aside block b
  datablock_test <- subset(datasubset, 
                           (yr>=blocks[[b]]$ymin) & 
                             (yr<=blocks[[b]]$ymax) )
  
  # keep the rest
  datablock_train <- subset(datasubset,
                            !((yr>=blocks[[b]]$ymin) & 
                                (yr<=blocks[[b]]$ymax) ))
  
  
  for (j in 1:n_zooClasses) {
    
    print(paste("zoo class:",zooClasses[j]))
    datablock_test$Z <- datablock_test[[zooClasses[j]]] # this could be log or not
    datablock_train$Z <- datablock_train[[zooClasses[j]]] # this could be log or not
    
    rf <- randomForest(Z ~ yday + sst + mo_sst_anom + sss + logchl + logdepth + SPG_proxy_anom_annual + logmld + bloom_start + bloom_duration + no3_winter + no3_mo_anom,
                       data=datablock_train,
                       importance=TRUE,
                       ntree=500,
                       mtry=opt_mtry_Z[j],
                       corr.bias=FALSE)
    print(rf)
    
    # get predictions and observations
    #    on training data (out of block; i.e. leave-one-out):
    pred_training <- predict(rf,datablock_train)
    obs_training <- datablock_train$Z
    #    on test data (in block):
    pred_test <- predict(rf,datablock_test)
    obs_test <- datablock_test$Z
    
    # get predictor importance
    varImpPlot(rf,
               type=1,
               scale=TRUE,
               cex=0.8,
               main=paste(zooClasses[j],"time block ",blocks[[b]]$ymin,"-",blocks[[b]]$ymax) )
    
    imp <-as.data.frame(importance(rf))
    imp <- imp[order(row.names(imp)), ]
    output_imp[[zooClasses[j]]][,(b+1)] <- imp$`%IncMSE`
    
    # plot observations vs predictions
    plot(obs_training,pred_training,main=paste("training: left out block",b,zooClasses[j]))
    plot(obs_test,pred_test,main=paste("test: block",b,zooClasses[j]))
    
    # get metrics... many metrics
    #get OOB error (mean square residuals), r^2 & variance explained
    output[[zooClasses[j]]]$OOB_error_training[b] <- rf$mse[length(rf$mse)]
    output[[zooClasses[j]]]$rOOB_error_training[b] <- sqrt(output[[zooClasses[j]]]$OOB_error_training[b])
    output[[zooClasses[j]]]$var_explained_training[b] <- rf$rsq[length(rf$rsq)]*100
    output[[zooClasses[j]]]$rmse_training[b]<- rmse(obs_training, pred_training) # are these equal?
    output[[zooClasses[j]]]$rsq_training[b] <- cor(obs_training, pred_training) ^ 2
    output[[zooClasses[j]]]$rmse_test[b]<- rmse(obs_test, pred_test)
    output[[zooClasses[j]]]$rsq_test[b] <- cor(obs_test, pred_test) ^ 2
  }
}

for (j in 1:n_zooClasses) {
  
  write.csv(output[[zooClasses[j]]],
            file=paste0("sensitivity analysis/","extrap_time_metrics_",zooClasses[j],".csv"),
            row.names = FALSE)
  write.csv(output_imp[[zooClasses[j]]],
            file=paste0("sensitivity analysis/","extrap_time_imp_",zooClasses[j],".csv"),
            row.names = FALSE)
  
  
}


