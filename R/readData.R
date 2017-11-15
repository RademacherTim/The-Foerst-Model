#============================================================================!
# Read all_sites_data.met and convert it to Rdata
#----------------------------------------------------------------------------!
library (R.matlab)
#library (tibble)

readDataSurv <- function (warning = 0) {
  data <- readMat ('matlab/all_sites_data_mixed.mat')
  dataSurv <<- data [["data.surv"]]
  return (dataSurv)
}

readDataDead <- function (warning = 0) {
  data <- readMat ('matlab/all_sites_data_mixed.mat')
  dataDead <<- data [["data.dead"]]
  return (dataDead)
}

#data$time
#============================================================================!
