validationNOTlongterm = 1; %1 = comparison to empirical data, 0 = long term simulation

harv_parameters = [0 1];
water_supplies = Inf;

load('all_sites_data_mixed');
start_age = 0; %age of the stand at the first observation (does not affect simulation); not known for mixed stand data

warning('off','all');

dx = 0.03; %resolution of stem size dimension, in meters
sizes = (dx:dx:1.5)'; %stem size scale, in metres
sizes2 = sizes.^2; %square of stem sizes
sizes3 = sizes.^3; %square of stem sizes
sizeclasses = length(sizes); %number of different stem sizes considered
dt = 0.1; %time step, in years

%the following is only for plotting
if validationNOTlongterm
    checkpoints = time-time(1); %points in simulation time at which empirical data is available
    checkpoints(4) = []; %ignore this one to allow 3x3 plot
    years = checkpoints(end); %number of years simulated
else        
    checkpoints = [time(end)-time(1);[125 200 300 500 1000]']; %last point for which empirical data is available, and several years after that
    years = checkpoints(end); %number of years simulated
end
timesteps = ceil(years/dt); %number of time steps in simulation

Light = zeros(sizeclasses,1);
TotalLeafArea = zeros(1,2);
GrowthFlux = zeros(sizeclasses,1);
GrowthGradient = zeros(sizeclasses,1);
DiameterGrowth = zeros(sizeclasses,1);
DeathRate  = zeros(sizeclasses,1);
u = zeros(sizeclasses,2);
u1 = zeros(sizeclasses,1);

%PARAMETERS
GrowthRate = [0.0055 0.0072]; %first entry: beech, second entry: spruce
maxsize = [2.7 1.7]; %from www.monumentaltrees.com
RespirationRate = GrowthRate./maxsize;
WaterCapacity = water_supplies;
lambda = [-0.042 -0.04]; %proportional to beer-lambert extinction coefficient
DeathPar = [0.000008 0.00005];
ReproductionRate = [5000 5000]; %estimated via Yoda's rule

HarvestIndex = ceil(harv_parameters(1)*length(sizes));
RegenYears = harv_parameters(2);
HarvestRate = 0;

%INITIAL CONDITIONS. Set initial u according to first empirical observation
for species=1:2
    u0 = 0.001*[data_surv{species}(:,1);data_dead{species}(:,1)]; %empirical diameter data is in mm, convert to meters
    u(1,species) = sum(u0 <= sizes(1));
    for i=2:length(sizes)
        k = (u0 > sizes(i-1)) & (u0 <= sizes(i));
        u(i,species) = round(sum(k)/dx); %convert number of individuals with diameter between sizes(i-1) and sizes (i) to density
    end
end

Fig1 = figure(1); hold on;
if validationNOTlongterm
    subplot(3,3,1); hold on;
    for species = 1:2
        v = zeros(length(sizes),1);
        v0 = 0.001*[data_surv{species}(:,1);data_dead{species}(:,1)];   

        v(1) = round(sum(v0 <= sizes(1))/dx);
        for i=2:length(sizes)
             k = (v0 > sizes(i-1)) & (v0 <= sizes(i));
             v(i) = round(sum(k)/dx); %initial empirical data (same as initial u)
        end
        if species == 1
            plot(sizes,v,':b','Linewidth',2);
            plot(sizes(6:end),u(6:end,species),'-b','Linewidth',1.5);
        else
            plot(sizes,v,':r','Linewidth',2);
            plot(sizes(6:end),u(6:end,species),'-r','Linewidth',1.5);
        end           
    end
    title('0 years');
    xlabel('dbh [m]');
end


for t=1:timesteps
   if ~ishghandle(Fig1)
        break;
   end
   
   %calculate portion of light (in [0,1]) arriving at heights corresponding
   %to different stem diameter sizes
   Light = ones(sizeclasses,1);
   for species=1:2
       Ind_Light = zeros(sizeclasses,1);
       for i=sizeclasses-1:-1:1
           Ind_Light(i) = Ind_Light(i+1) + u(i+1,species) * sizes2(i+1);
       end
       TotalLeafArea(species) = dx*(Ind_Light(1) + u(1)*sizes2(1));
       Ind_Light = Ind_Light + 0.5 * u(:,species) .* sizes2;

       Light = Light .* exp(lambda(species) * dx * Ind_Light);
   end

   %gross photosynthetic production assuming no water constraints
   %if this exceeds WaterCapacity, then all actual growth is reduced
   %accordingly
   PotentialProduction = 0;
   for species=1:2
       PotentialProduction = PotentialProduction + sum(GrowthRate(species) * Light .* u(i,species) .* sizes2) * dx;
   end

   for species=1:2 %1 = beech, 2 = spruce
       
       %next, we calculate the different components determining the size
       %distribution one time step later. this is done as described in the
       %paper
       
       %Growth
       DiameterGrowth = (GrowthRate(species) * Light * min(1,WaterCapacity/PotentialProduction) - RespirationRate(species) * sizes);

       %Death
       MassGrowth = DiameterGrowth .* sizes2;
       DeathRate = 0.5*DeathPar(species)./max(0,MassGrowth);
       DeathRate(isnan(DeathRate)) = 0;

       GrowthFlux = DiameterGrowth .* u(:,species);
       GrowthGradient(1) = (GrowthFlux(2)-GrowthFlux(1))/dx;
       GrowthGradient(end) = (GrowthFlux(end)-GrowthFlux(end-1))/dx;
       GrowthGradient(2:end-1) = (GrowthFlux(3:end)-GrowthFlux(1:end-2))/(2*dx);

       %PDE
       u1(:) = max(0,u(:,species) - dt * (GrowthGradient + (DeathRate + HarvestRate).* u(:,species))); %"max(0," isnt really necessary here, but an efficient way to force stability
       u1(1) = ReproductionRate(species) * TotalLeafArea(species); %birth of new trees
       u(:,species) = u1;
   end

   %PLOTTING
   k0 = find(t*dt == checkpoints); %checkpoints[k0] = t*dt, if there is such a k0
   if ~isempty(k0)
       if validationNOTlongterm
           subplot(3,3,k0); hold on;

           for species = 1:2
               v = zeros(length(sizes),1);                
               v0 = 0.001*[data_surv{species}(:,k0);data_dead{species}(:,k0)]; %empirical data

               v(1) = round(sum(v0 <= sizes(1))/dx);
               for i=2:length(sizes)
                    k = (v0 > sizes(i-1)) & (v0 <= sizes(i));
                    v(i) = round(sum(k)/dx); %convert absolute individuals to density
               end
               if species == 1
                   plot(sizes,v,':b','Linewidth',2);
                   plot(sizes(6:end),u(6:end,species),'-b','Linewidth',1.5);
               else
                   plot(sizes,v,':r','Linewidth',2);
                   plot(sizes(6:end),u(6:end,species),'-r','Linewidth',1.5);
               end           
           end
           title(strcat(num2str(start_age+round(t*dt)),' years'));
           xlabel('dbh [m]');
       else
           figure(1); hold on;
           plot3(checkpoints(k0)*ones(length(sizes(6:end)),1),sizes(6:end),u(6:end,1),'-b','Linewidth',1.5); %first dimension is simply time, to have a nice 3D plot
           plot3(checkpoints(k0)*ones(length(sizes(6:end)),1),sizes(6:end),u(6:end,2),'-r','Linewidth',1.5);
           set(gca,'xscale','log');
           ylabel('dbh [m]');
           xlabel('age [y]');
           xlim([checkpoints(1) checkpoints(end)]); 
       end
   end
end
 
 
 