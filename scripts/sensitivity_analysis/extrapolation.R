# apply block cross validation in space to test ability to extrapolate
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

# --- design the blocks
# plot samples by species
plot(dataset$x,dataset$y,main="CPR samples and C fin>0",col="grey",xlab="lon",ylab="lat")
points(subset( dataset, Cfin>0, select=c(x,y)),col="black")

plot(dataset$x,dataset$y,main="CPR samples and C hel>0",col="grey",xlab="lon",ylab="lat")
points(subset( dataset, Chel>0, select=c(x,y)),col="black")

plot(dataset$x,dataset$y,main="CPR samples and small copepods>0",col="grey",xlab="lon",ylab="lat")
points(subset( dataset, small>0, select=c(x,y)),col="black")  

# define models
mdl_name <- c("NEAtlantic") 
# subsample for full model
doSubset <-  FALSE
subset_size <- 0.5
if (doSubset & mdl_name=="NEAtlantic") {
  subsample_rows <- sample(1:n_samples,subset_size*n_samples) # get random sub-sample of total obs
  datasubset <- dataset[subsample_rows,] 
} else {
  datasubset <- dataset
}

print(paste("datasubset has",nrow(datasubset),"records"))

rf_string <- "full_with_no3"
datasubset <- datasubset[datasubset$sss>28,] # take out low salinities 

# --- define spatial blocks - coarsely to start with
# (can reference using blocks$block1$E etc)
blocks <- list(
  block1=data.frame(W =-20, E = -10, S=40,N=48),
  block2=data.frame(W =-20, E = -10, S=48,N=56),
  block3=data.frame(W =-20, E = -10, S=56,N=65),
  block4=data.frame(W =-10, E = 0, S=40,N=48),
  block5=data.frame(W =-10, E = 0, S=48,N=56),
  block6=data.frame(W =-10, E = 0, S=56,N=65),
  block7=data.frame(W =0, E = 10, S=48,N=56),
  block8=data.frame(W =0, E = 10, S=56,N=65)
)

n_blocks <- length(blocks)
blocks_col <-rainbow(n_blocks) # define plotting colours

#plot blocks
plot(datasubset$x,datasubset$y)
for (b in 1:n_blocks) {
  datablock_in <- subset(datasubset, 
                         (x>blocks[[b]]$W) & 
                           (x<blocks[[b]]$E) & 
                           (y>blocks[[b]]$S) & 
                           (y<blocks[[b]]$N)  )
  print(paste('block',b,'has',nrow(datablock_in),'samples'))

  points(datablock_in$x,datablock_in$y,col=blocks_col[b])
}

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

# --- loop over blocks
for (b in 1:n_blocks) {
  
  # put aside block b
  datablock_test <- subset(datasubset, 
                        (x>blocks[[b]]$W) & 
                        (x<blocks[[b]]$E) & 
                        (y>blocks[[b]]$S) & 
                        (y<blocks[[b]]$N)  )
  
  # keep the rest
  datablock_train <- subset(datasubset,
                          !((x>blocks[[b]]$W) & 
                            (x<blocks[[b]]$E) & 
                            (y>blocks[[b]]$S) & 
                            (y<blocks[[b]]$N)))


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

# now get RMSE and R^2 for full data
rmse_full <- rep(NA, 3)
rsq_full <- rep(NA, 3)
for (j in 1:n_zooClasses) {
  
  print(paste("zoo class:",zooClasses[j]))
  datasubset$Z <- datasubset[[zooClasses[j]]]
  rf <- randomForest(Z ~ yday + sst + mo_sst_anom + sss + logchl + logdepth + SPG_proxy_anom_annual + logmld + bloom_start + bloom_duration + no3_winter + no3_mo_anom,
                     data=datasubset,
                     importance=TRUE,
                     ntree=500,
                     mtry=opt_mtry_Z[j],
                     corr.bias=FALSE)
  print(rf)
  
  # get predictions and observations
  pred_training <- predict(rf,datasubset)
  obs_training <- datasubset$Z
  
  # get rmse and rsq
  rmse_full[j] <- rmse(obs_training, pred_training) # are these equal?
  rsq_full[j] <- cor(obs_training, pred_training) ^ 2

}





for (j in 1:n_zooClasses) {
  
  write.csv(output[[zooClasses[j]]],
          file=paste0("sensitivity analysis/","extrap_metrics_",zooClasses[j],".csv"),
          row.names = FALSE)
  
}


  