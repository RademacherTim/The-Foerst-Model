#============================================================================!
# Read empirical data for beech and spruce from the permanent sampling plots 
# surveyed by the Technische Universitaät München (TUM) from 1848? to present.
# 
# The data is used for evaluation of the först model (först_model.R) of forest 
# size distributions.
#----------------------------------------------------------------------------!
library (R.matlab)
#library (tibble)

readDataSurv <- function (warning = 0) {
  data <- readMat ('matlab/all_sites_data_mixed.mat')
  dataSurv <- data [["data.surv"]]
  time <- data [['time']]
  return (dataSurv)
}

readDataTime <- function (warning = 0) {
  data <- readMat ('matlab/all_sites_data_mixed.mat')
  dataTime <- data [["time"]]
  return (dataTime)
}

readDataDead <- function (warning = 0) {
  data <- readMat ('matlab/all_sites_data_mixed.mat')
  dataDead <- data [["data.dead"]]
  return (dataDead)
}

#data$time
#============================================================================!
