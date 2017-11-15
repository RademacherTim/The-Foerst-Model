#============================================================================!
# 
#----------------------------------------------------------------------------!

# start from clean slate
#----------------------------------------------------------------------------!
rm (list = ls ())

# Load necessary dependencies
#----------------------------------------------------------------------------!
source ('R/readData.R')

water_supplies   <- Inf      #
harv_parameters <- c (0, 1) # 
nspecies <- 2               # Beech is species = 1, spruce is species = 2

start_age   <- 0      # age of the stand at the first observation (does not 
                      # affect simulation); not known for mixed stand data
  
nspecies    <- 2                      # number of species in simulation
dx          <- 0.03                   # resolution of stem size dimension [m]
sizes       <- seq (dx, 1.5, by = dx) # stem size scale, in metres
sizes2      <- sizes^2                # square of stem sizes
sizes3      <- sizes^3                # cube of stem sizes
sizeclasses <- length (sizes)         # number of different stem sizes considered
dt          <- 0.01                   # time steps per year
years       <- 1000                   # number of years in simulation
timesteps   <- ceiling (years / dt)   # number of time steps in simulation

# Set initial values to zero
#----------------------------------------------------------------------------!
Light          <- rep (0, sizeclasses)
TotalLeafArea  <- rep (0, nspecies)
GrowthFlux     <- rep (0, sizeclasses)
GrowthGradient <- rep (0, sizeclasses)
DiameterGrowth <- rep (0, sizeclasses)
DeathRate      <- rep (0, sizeclasses)
u              <- matrix (rep (0, nspecies * sizeclasses), nrow = sizeclasses)
u1             <- rep (0, sizeclasses)

# Species parameters for 
#     - beech  (first  entry)
#     - spruce (second entry)
#----------------------------------------------------------------------------!
BeechParameters  <- c (0.0055, 2.7, -0.042, 0.000008, 5000.0)
SpruceParameters <- c (0.0072, 1.7, -0.040, 0.000050, 5000.0)
GrowthRate       <- c (BeechParameters [1], SpruceParameters [1]) # species-specific growth rate []
maxsize          <- c (BeechParameters [2], SpruceParameters [2]) # max tree size from www.monumentaltrees.com
RespirationRate  <- GrowthRate / maxsize                          # respiration rate []
WaterCapacity    <- water_supplies                                #
lambda           <- c (BeechParameters [3], SpruceParameters [3]) # proportional to beer-lambert extinction coefficient
DeathPar         <- c (BeechParameters [4], SpruceParameters [4]) #
ReproductionRate <- c (BeechParameters [5], SpruceParameters [5]) # estimated via Yoda's rule

# harvest parameters
#----------------------------------------------------------------------------!
HarvestIndex <- ceiling (harv_parameters [1] * sizeclasses)
RegenYears   <- harv_parameters [2]
HarvestRate  <- 0.0

# set initial size distribution (u) according to first observation
#----------------------------------------------------------------------------!
dataSurv <- readDataSurv ()
dataDead <- readDataDead ()
for (species in 1:nspecies) {
  # empirical diameter data is in mm, convert to meters
  u0 <- 0.001 * c (dataSurv [[species]] [[1]] [, 1], dataDead [[species]] [[1]] [, 1])
  u [1, species] <- sum (u0 <= sizes [1])
  for (i in 2:sizeclasses) {
    k <- (u0 > sizes [i - 1]) & (u0 <= sizes [i])
    u [i, species] <- round (sum (k) / dx) # convert number of individuals with diameter between sizes(i-1) and sizes (i) to density
  }
}

# Set timer to measure execution time
#----------------------------------------------------------------------------!
start.time <- Sys.time ()

# Loop over timesteps
#----------------------------------------------------------------------------!
for (t in 1:timesteps) {
  # calculate portion of light (Light) arriving at heights corresponding
  # to different stem diameter sizes
  #--------------------------------------------------------------------------!
  Light <- rep (1, sizeclasses) # Initially assume all are fully illuminated
  for (species in 1:nspecies) {
    Ind_Light <- rep (0, sizeclasses)
    for (i in (sizeclasses-1):1) {
      Ind_Light [i] <- Ind_Light [i + 1] + u [i + 1, species] * sizes2 [i + 1]
    }
    TotalLeafArea [species] <- dx * (Ind_Light [1] + u [1] * sizes2 [1])
    Ind_Light <- Ind_Light + 0.5 * u [, species] * sizes2
    Light <- Light * exp (lambda [species] * dx * Ind_Light)
  }

  # gross photosynthetic production assuming no water constraints, if this 
  # exceeds WaterCapacity, then all actual growth is reduced accordingly
  PotentialProduction <- 0
  for (species in 1:nspecies) {
    PotentialProduction <- PotentialProduction + sum (GrowthRate [species] * Light * u [i, species] * sizes2) * dx  
  }

  for (species in 1:nspecies) {

    # calculate the different components determining the size distribution
    # one time step later. this is done as described in the paper

    # Growth
    DiameterGrowth = (GrowthRate [species] * Light * min (1, WaterCapacity / 
                                                             PotentialProduction) 
                     - RespirationRate [species] * sizes)

    # Death
    MassGrowth <- DiameterGrowth * sizes2
    DeathRate  <- 0.5 * DeathPar [species] / max (0, MassGrowth)
    DeathRate [is.nan (DeathRate)] <- 0

    end <- sizeclasses
    GrowthFlux                 <- DiameterGrowth * u [, species]
    GrowthGradient [1]         <- (GrowthFlux [2]     - GrowthFlux [1])         / dx
    GrowthGradient [end]       <- (GrowthFlux [end]   - GrowthFlux [end-1])     / dx
    GrowthGradient [2:(end-1)] <- (GrowthFlux [3:end] - GrowthFlux [1:(end-2)]) / (2 * dx);

    #PDE
    #"max(0," isnt really necessary here, but an efficient way to force stability
    u1     <- max (0, u[ , species] - dt * (GrowthGradient + (DeathRate + HarvestRate) * u [, species])) 
    u1 [1] <- ReproductionRate [species] * TotalLeafArea [species] # birth of new trees
    u [, species] <- u1
  }
}

# Determine time it took to run
#----------------------------------------------------------------------------!
end.time   <- Sys.time ()
time.taken <- end.time - start.time
time.taken
#============================================================================!