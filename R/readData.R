#============================================================================!
# Read all_sites_data.met and convert it to Rdata
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
