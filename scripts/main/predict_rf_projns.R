# --- Take fitted RF model(s) and predict on projection time steps ---
# -----------------------------------------
# EJT December 2025 - tidied up for manuscript submission

# Note: predicted on monthly mean bias-corrected regional ocean projections.
# Original data (not bias corrected) are available from BODC at
# https://doi.org/10.5285/07877700-0e22-5fb5-e063-6c86abc03058 and
# https://doi.org/10.5285/0786d770-ba54-0e57-e063-6c86abc09fdf, and through
# Zendo at https://zenodo.org/records/3953801.

# Please contact the author for bias corrected datafiles.

# required packages - do I need all these?
library(randomForest)
library(ggplot2)
library(dplyr) 
library(lubridate)

rm(list = ls())

# ---- load cleaned ZED and AMM7 dataset ---
# zed taxa are all log transformed using offset equal to half min non-zero value (taxon specific)
load("data/logzed_and_predictors.RData")
print("loaded logzed dataset")

# ---- define zooplankton classes
zooClasses <- c("small","Cfin","Chel")
n_zooClasses <- length(zooClasses)

# ---- Define model:
model <- c("NEAtlantic")  
projn <- c("RECICLE_GFDL") #,"RECICLE_IPSL","RECICLE_GFDL","HADGEM" # "reanalysis" # "HADGEM" # ocean model

# note that reduced/full_with_no3 is latest. need to use reduced for HADGEM because no bloom phenology.
if (projn=="HADGEM") { # random forest model
  rf_string <- "reduced_with_no3" #'final_with_yday'
} else {
  rf_string <- "full_with_no3"
}
print(paste("Using **", rf_string,"** RF model trained on data for **",model,"**"))

# region to predict on (NB can be same as model)
region <- c("NEAtlantic")  # predict on whole region
print(paste("Predicting on: **",region,"**"))

# ---- correct for bias?
correct_bias <- TRUE
if (correct_bias) {
  ROE_values <- read.csv(paste0("output/bias correction coefficients/ROE_values",model,"_",rf_string,".csv"))
  print(paste("loaded bias correction values for model: ",model))
}

# --- project on AMM7 projections using list of files ---

# monthly bias corrected projns:
# where to put these?
datapath <- "C:/Users/xqb21120/OneDrive - University of Strathclyde/Documents/MATLAB/ECOWINGS/randomForest/bias_corrected_monthly_projns/"

# where to store results:
outpath <- "output/projections/"

# time period
dt <- 'monthly' # 'daily'
decade<-c("2010s","2040s","2090s")
n_decades<-length(decade)

# pre load the RF models - - access using rf_all[[j]]
rf_by_zoo <- vector(mode = "list", length = 18)
rf_all <- list(rf_by_zoo,rf_by_zoo,rf_by_zoo)
for (j in 1:n_zooClasses) { # small, Cfin, Chel
  print('loading rf model for each zooplankton class...')
  print(zooClasses[j])
  mdl_name <- paste0('output/rf_models/rf_',rf_string,'_',model,'_',zooClasses[j],'.Rdata') # rf model to load
  print(mdl_name)
  load(mdl_name)
  rf_all[[j]] <- rf
}

# loop over time periods (now, near future, far future)
for (i in 1:n_decades) { # 1:n_decades
  
  print(paste("getting datafiles for ",decade[i]))
  searchStr <- glob2rx(paste0('*',projn,'_mo_NEAtlantic_',substring(decade[i],1,3),'*')) # make sure to exclude partial
  
  filenames <- list.files(datapath, pattern=searchStr)
  n_files <- length(filenames)
  
  # loop over time steps
  for (k in 1) { #:n_files) {
    
    fname <- paste0(datapath,filenames[k])
    print(filenames[k])
    
    # get datestamp from filename
    dateStr <- substr(fname,nchar(fname)-9,nchar(fname)-4)
    
    # read ocean model data for timestep
    data_timestep <- read.csv(fname, header = TRUE)
    
    # give variables the right names
    data_timestep <- data_timestep %>% rename(SPG_proxy_anom_annual = SPG_proxy_annual)
    data_timestep$logdepth <- log10(data_timestep$depth)
    
    idx<-which(data_timestep$mld<0)
    data_timestep$mld[idx]<-1 
    data_timestep$logmld <- log10(data_timestep$mld)
    
    # log transform chl: log10(chl+offset), where offset=chl_add=0.5*min(chl>0), as calculated by clean_up_data.R
    data_timestep$logchl <- log10(data_timestep$chl + chl_add)
    
    # --make empty results dataframe--
    ZED_predictions <- data.frame(t= data_timestep$t,
                                  x=data_timestep$x,
                                  y= data_timestep$y,
                                  depth=data_timestep$depth,
                                  boxmask=data_timestep$boxmask) # will make it easier to process the results!!
    n<-nrow(ZED_predictions)
    blank <- double(n)
    for (j in 1:n_zooClasses) {
      ZED_predictions[[zooClasses[j]]]<-blank
    }
    
    # --- loop over Z classes ---
    for (j in 1:n_zooClasses) {
      
      # get projections for this zoo class
      projected_ZE <- predict(rf_all[[j]], data_timestep)
      
      # bias correct
      if (correct_bias) {
        m <- ROE_values$m_ROE[ROE_values$zooClasses==zooClasses[j]]
        b <- ROE_values$b_ROE[ROE_values$zooClasses==zooClasses[j]]
        projected_ZE <- m*projected_ZE+b
      }
      
      # add to dataframe
      ZED_predictions[,zooClasses[j]] <- projected_ZE 
      
    }
    
    #---  save as csv ---
    fname = paste0(outpath,'predicted_ZE_',rf_string,'_',model,'_',region,'_',projn,'_',dt,'_',dateStr,'.txt')
    write.csv( format(ZED_predictions),
               file = fname,
               row.names=FALSE,
               quote=FALSE)
  }
}

# ---------- all done :) -----------