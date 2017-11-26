!============================================================================!
subroutine forst_model (u0,       & ! initial size distribution
                        ut,       & ! final size distribution
                        years,    & ! number of years in simulation
                        !dt,       & ! time steps per year
                        nspecies  & ! number of species in simulation

                        )
!============================================================================!
! The först model simulates forest size distribution employing a set of 
! partial differential equations (PDE)
!----------------------------------------------------------------------------!
!     developed by:
!     Robert Beyer    (rb792@cam.ac.uk)
!     Tim Rademacher  (trademacher@fas.harvard.edu)
!    
!     initially parameterised and evaluated against permanent plot 
!     measurements of beech and spruce from the Technische Universität 
!     München (TUM) from 1848? to present.
!----------------------------------------------------------------------------!

! start from clean slate
!----------------------------------------------------------------------------!
implicit none

! declare and initialise variables
!----------------------------------------------------------------------------!
integer :: i                                  ! loop variable
integer :: species                            ! loop variable
integer, intent (in) :: nspecies              ! number of species in simulation
real (kind = 8), parameter :: dx      = 0.03  ! resolution of stem size dimension [m]
integer, parameter :: nsizes = int (dx / 1.5) ! number of size classes up to 1.5m     
real (kind = 8), dimension (nsizes) :: sizes  = (/ (real (i*dx), i = 1, nsizes) /) ! Different size classes
real (kind = 8), dimension (nsizes) :: sizes2 = (/ (sizes (i)**2.0, i = 1, nsizes) /) ! Size classes squared
real (kind = 8), dimension (nsizes) :: sizes3 = (/ (sizes (i)**3.0, i = 1, nsizes) /) ! Size classes cubed
real (kind = 8), dimension (nspecies, nsizes) :: u0, u1, u2, ut ! initial, transient and final size distributions
!----------------------------------------------------------------------------!
integer, intent (in) :: years                ! number of years in simulation
real (kind = 8), parameter :: dt = 0.01      ! time step as proportion of year
integer :: t                                 ! loop variable for the time step 
integer :: timesteps                         ! number of time steps in simulation

real (kind = 8), dimension (nsizes) :: light ! proportion of light penetrating to each size class
real (kind = 8), dimension (nsizes) :: indLight !

real (kind = 8), dimension (5), parameter :: beechParameters  = (/0.0055, 2.7, -0.042, 0.000008, 5000.0/)
real (kind = 8), dimension (5), parameter :: spruceParameters = (/0.0072, 1.7, -0.040, 0.000050, 5000.0/)
real (kind = 8), dimension (nspecies) :: totalLeafArea      !
real (kind = 8), dimension (nspecies) :: lambda     		!
real (kind = 8), dimension (nspecies) :: growthRate         ! 
real (kind = 8), dimension (nspecies) :: maxSize            ! maximum tree stem diameter [m]
real (kind = 8), dimension (nspecies) :: respirationRate    !
real (kind = 8), dimension (nspecies) :: deathPar           !
real (kind = 8), dimension (nspecies) :: reproductionRate   !

real (kind = 8), parameter :: waterSupplies = huge (0.0d0) ! Set water supply to infinity
real (kind = 8), dimension (nsizes) :: waterCapacity   = waterSupplies
real (kind = 8), dimension (nsizes) :: growthFlux      = 0 !
real (kind = 8), dimension (nsizes) :: growthGradient  = 0 !
real (kind = 8), dimension (nsizes) :: diameterGrowth  = 0 !
real (kind = 8), dimension (nsizes) :: deathRate       = 0 !
real (kind = 8), dimension (nsizes) :: massGrowth      = 0 !
real (kind = 8) :: potentialProduction                !

! harvest parameters
!----------------------------------------------------------------------------!
real (kind = 8), dimension (2), parameter :: harvParameters = (/0, 1/)
real (kind = 8), parameter :: harvestIndex = ceiling (harvParameters (1) * nsizes) ! nsizes or size classes
real (kind = 8), parameter :: regenYears   = harvParameters (2)
real (kind = 8), parameter :: harvestRate = 0.0

! initialise variables
!----------------------------------------------------------------------------!
totalLeafArea = 0

! set species-specific parameters
!----------------------------------------------------------------------------!
if (nspecies == 1) then
  growthRate       = beechParameters (1)
  maxSize          = beechParameters (2)
  lambda           = beechParameters (3)
  deathPar         = beechParameters (4)
  reproductionRate = beechParameters (5)
else if (nspecies == 2) then
  growthRate       = (/beechParameters (1), spruceParameters (1)/)
  maxSize          = (/beechParameters (2), spruceParameters (2)/)
  lambda           = (/beechParameters (3), spruceParameters (3)/)
  deathPar         = (/beechParameters (4), spruceParameters (4)/)
  reproductionRate = (/beechParameters (5), spruceParameters (5)/)
end if
respirationRate = growthRate / maxSize

! calculate the number of time steps
!----------------------------------------------------------------------------!
timesteps = ceiling (years / dt)

! set initial size distribution
!----------------------------------------------------------------------------!
u1 = u0

! loop over timesteps
!----------------------------------------------------------------------------!
time : do t = 1, timesteps
  ! calculate portion of light (Light) arriving at heights corresponding
  ! to different stem diameter sizes
  !--------------------------------------------------------------------------!
  light = 1    ! Initially assume all are fully illuminated
  calcLight : do species = 1, nspecies
    indLight = 0 ! Set light penetration to 0  
    layers : do  i = nsizes-1, 1, -1 ! Loop over canopy layers
      indLight (i) = indLight (i + 1) + u1 (i + 1, species) * sizes2 (i + 1)
    end do layers ! loop ending
    totalLeafArea (species) = dx * (indLight (1) + u1 (1, species) * sizes2 (1)) ! Not sure u is correct here 
    indLight = indLight + 0.5 * u1 (:, species) * sizes2
    light    = light * exp (lambda (species) * dx * indLight)
  end do calcLight ! loop ending

  ! gross photosynthetic production assuming no water constraints, if this 
  ! exceeds WaterCapacity, then all actual growth is reduced accordingly
  !--------------------------------------------------------------------------!
  potentialProduction = 0
  calcPotentialProduction : do species = 1, nspecies
    potentialProduction = potentialProduction +                              &
                          sum (growthRate (species) * light *                &
                               u1 (:, species) * sizes2) * dx  
  end do calcPotentialProduction ! loop ending

  calcSizeDistribution : do species = 1, nspecies
    ! calculate the different components determining the size distribution
    ! one time step later. this is done as described in the paper

    ! Growth
    !------------------------------------------------------------------------!
    diameterGrowth = (growthRate (species) * light *                         &
                      min (1.0, waterCapacity / potentialProduction) -       &
                      respirationRate (species) * sizes)

    ! Death
    massGrowth = diameterGrowth * sizes2
    deathRate  = 0.5 * deathPar (species) / max (0.0, massGrowth)
    where (isnan (deathRate)) deathRate = 0

    growthFlux                    = diameterGrowth * u1 (:, species)
    growthGradient (1)            = (growthFlux (2)        - growthFlux (1))         / dx
    growthGradient (nsizes)       = (growthFlux (nsizes)   - growthFlux (nsizes-1))     / dx
    growthGradient (2:(nsizes-1)) = (growthFlux (3:nsizes) - growthFlux (1:(nsizes-2))) / (2 * dx);

    !Partial Differential Equation
    !"max(0," isnt really necessary here, but an efficient way to force stability
    u2              = max (0, u (:, species) - dt * (GrowthGradient + (DeathRate + HarvestRate) * u1 (:, species))) 
    u2 (1, species) = ReproductionRate (species) * TotalLeafArea (species) ! birth of new trees
    u1              = u2
  end do calcSizeDistribution ! loop ending

  ut = u1
  write (*, *) 'Done with loop ',t
end do time ! loop ending

!============================================================================!
end subroutine
!============================================================================!