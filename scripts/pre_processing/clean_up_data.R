# --- Load CPR samples and predictors from amm7 and climate indices ---
# Clean data ready for RF modelling.
# -----------------------------------------
# EJT December 2025 - tidied up for manuscript submission

# Note: dataset "zed_with_amm7" is derived by matching CPR samples with the
# phys-bgc ocean reanalysis data. CPR data available from the Marine Biological
# Association at https://doi.dassh.ac.uk/data/1850 and
# https://doi.mba.ac.uk/data/3118. physical-biogeochemical model reanalysis data
# available through Copernicus Marine Service
# (https://data.marine.copernicus.eu), products NWSHELF_MULTIYEAR_PHY_004_009
# and NWSHELF_MULTIYEAR_BGC_004_011.

# required packages
library(randomForest)
library(ggplot2)
library(dplyr) 
library(lubridate)

# ---- load data ----
zed_with_amm7 <- read.csv("C:/Users/xqb21120/OneDrive - University of Strathclyde/Documents/MATLAB/ECOWINGS/randomForest/zed_and_amm7.txt", header = TRUE)
zed_with_amm7 <- zed_with_amm7[zed_with_amm7$yr>=1998,] # chuck out pre 1998 when chl not useful

# additional predictor: SPG proxy based on salinity time series having strongest correlation with SPG index
SPG_proxy <- read.csv("data/SPG_proxy_annual.txt")

# merge with amm7
zed_with_amm7$month <- month(dmy(zed_with_amm7$t))

zed_with_amm7 <- merge(zed_with_amm7,SPG_proxy,by="yr",all.x=TRUE)

# rf can't handle nans so chuck these rows out. corresponds to samples outwith amm7 domain.
# only 53075 samples out of 246611 are within amm7 domain = 22%
zed_with_amm7 <- na.omit(zed_with_amm7)

# check variables correctly classified (e.g. categorical as factor, character as date)
zed_with_amm7$is_day <- as.factor(zed_with_amm7$is_day)

# remove cpr chl and rename amm7_chl
zed_with_amm7 <- subset(zed_with_amm7, select = -c(chl) )
zed_with_amm7 <- zed_with_amm7 %>% rename(chl = chl_amm7)

# log transform chl
chl_add <- 0.5*min(zed_with_amm7$chl[zed_with_amm7$chl>0]) # nominal value to add to chl
zed_with_amm7$logchl <- log10( (zed_with_amm7$chl) + chl_add )

# ---- set zooplankton groups ----
# groups to model are: oithona, appendicularia, small copepods, psuedocalanus, C fin, C hel, large copepods, euphausiids

# ------ define "small copepod prey"
# initially tried two separate definitions, based on Olin et al. (2022). https://doi.org/10.1093/icesjms/fsac101 & MacDonald, et al. (2018). https://doi.org/10.3389/fmars.2018.00339
# didn't affect results. went with MacDonald definition.

# MacDonald: C unid I-II, Centropages hamatus, Centropages typicus, Pseudocalanus, Temora longicornis, Acartia spp., Parapseudocalanus, Oithona, Microcalanus
# Olin: Acartia spp., Oithona spp., Para-Pseudocalanus spp., Temora longicornis

# Olin defn:
zed_with_amm7 <- zed_with_amm7 %>% rename(small_Olin = small)

# MacDonald defn:
zed_with_amm7$small <- rowSums(zed_with_amm7[ ,
                                                        c("C14",
                                                          "centro_hammatus",
                                                          "centro",
                                                          "centro_typicus",
                                                          "parapseudo",
                                                          "temora",
                                                          "acartia",
                                                          "oith")])


# large copepods = sum of "Cfin" and "Chel"
zed_with_amm7$large_copepods <- rowSums(zed_with_amm7[ ,
                                                       c("Cfin",
                                                         "Chel",
                                                         "C56unid")])

# "other food"
zed_with_amm7$other <- rowSums(zed_with_amm7[ ,c("appen",
                                                 "cirri",
                                                 "nauplii",
                                                 "decapoda",
                                                 "evadne",
                                                 "eggs",
                                                 "fish",
                                                 "hyperiidea",
                                                 "metridia",
                                                 "podon")])   

                                                 

# get log depth & mld
zed_with_amm7$logdepth <- log10(zed_with_amm7$depth)
zed_with_amm7$logmld <- log10(zed_with_amm7$mld)

# log transform all taxa, using half min as offset
# this includes "total"
taxa <- names(zed_with_amm7)[c(6:30,56:58)]
logzed_offsets <- double(length(taxa))
logzed_with_amm7 <- zed_with_amm7
for (i in 1:length(taxa)) {
  t<-taxa[i]
  print(t)
  idx <- zed_with_amm7[[t]]>0
  Z_add <- 0.5*min( zed_with_amm7[[t]][idx] ) # add half min of non zero values
  logzed_offsets[i] <- Z_add
  print(Z_add)
  logzed_with_amm7[[t]] <- log10( zed_with_amm7[[t]] + Z_add )
}

write.csv(data.frame(taxa,logzed_offsets),file="data/logzed_offsets.csv",row.names = FALSE)

# ---- save for reading into final models stage ---- 
save(zed_with_amm7, file = "data/zed_and_predictors.RData")
save(logzed_with_amm7, logzed_offsets, chl_add,file = "data/logzed_and_predictors.RData")

# also save in csv file in case Matlab wants it
write.csv(logzed_with_amm7,file="data/logzed_and_predictors.csv",row.names = FALSE)
write.csv(zed_with_amm7,file="data/zed_and_predictors.csv",row.names = FALSE)
