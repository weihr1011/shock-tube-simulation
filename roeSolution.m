%roeSolution.m
%Purpose: uses Godunov's scheme and Roe's algorithm to numerically   
%         simulate the time evolution of Sod's shock tube 
%Parameters: (passed through par structure)
%            densL       initial density in the left chamber
%            vxL         initial velocity in the left chamber
%            presL       initial pressure in the left chamber
%            densR       initial density in the right chamber
%            vxR         initial velocity in the right chamber
%            presR       initial pressure in the right chamber
%            x0          the diaphragm position
%            cellCoords  the x-positions of the cells in the 
%                        shock tube grid
%            gamma       the specific heat ration for the ideal gas
%            t           the time at which the function will terminate
%                        the simulation and return the new state of
%                        the system
%            maxCycles   the maximum number of times that the simulation 
%                        will apply Godunov's schem (prevents infinite
%                        loops if the timesteps become too small to
%                        reach tLim)
%            cfl         the Courant-Friedrichs-Lewy number for 
%                        determining timestep size
%Outputs: density  a vector containing the final result for the gas
%                  density in every cell (same size as cellCoords param)
%         velocity a vector containing the final result for the gas
%                  velocity in every cell (same size as cellCoords param) 
%         pressure a vector containing the final result for the gas
%                  pressure in every cell (same size as cellCoords param)
%


function [density,velocity,pressure] = roeSolution(par)

%Before starting the simulation, we need to create a grid/vector that
%matches the initial state of the shock tube system. 

%use the list of x-positions to determine the size
%of the simulation grid
numXCells = size(par.cellCoords, 2);
%only run the simulation if there are at least two cells
assert(numXCells > 1);
%calculate the cell width (could also just use par.dx parameter)
dx = par.cellCoords(1,2) - par.cellCoords(1,1);

%determine the index for the boundary between
%the left and right chambers
diaphragmIndex = ceil((par.x0-par.cellCoords(1,1))/dx);

%create the density vector
leftDens = par.densL*ones([1 diaphragmIndex]);
rightDens = par.densR*ones([1 numXCells-diaphragmIndex]);
currDens = [leftDens rightDens];

%create the velocity vector
leftVel = par.vxL*ones([1 diaphragmIndex]);
rightVel = par.vxR*ones([1 numXCells-diaphragmIndex]);
velocity = [leftVel rightVel];
%create the pressure vector
leftPres = par.presL*ones([1 diaphragmIndex]);
rightPres = par.presR*ones([1 numXCells-diaphragmIndex]);
pressure = [leftPres rightPres];
%create a 3xnumXCells matrix that stores the conservative variables
%for each x-position (density, momentum, and energy, in that order)
momentum = velocity.*currDens;
energy = 0.5*currDens.*velocity.^2+pressure./(par.gamma-1);
U(1,:) = currDens;
U(2,:) = momentum;
U(3,:) = energy;

%start the simulation at t=0
t = 0.0;

%keep track of how many times we have updated the simulation grid
cycle = 0;

%run the simulation until we reach the desired time 
%or the maximum number of cycles
while t < par.t && cycle <= par.maxCycles
    
    %We first add a ghost cell to each side of the state vectors.
    %Since we use transparent boundary conditions, we simply have to 
    %copy the state variables in the left and rightmost cells.
    
    %create density vector with ghost cells
    currDens = [currDens(1,1) currDens currDens(1,size(currDens,2))];
    %copy the densities in every i-th cell 
    %(i.e., the left cell densities)
    dens_i = currDens(1,1:numXCells+1);
    %copy the densities in every i+1-th cell 
    %(i.e., the right cell densities)
    dens_i1 = currDens(1,2:numXCells+2);
    
    %create velocity vector with ghost cells
    velocity = [velocity(1,1) velocity velocity(1,size(velocity,2))];
    %copy the velocities in every i-th cell 
    %(i.e., the left cell velocities)
    vel_i = velocity(1,1:numXCells+1);
    %copy the velocities in every i+1-th cell 
    %(i.e., the right cell velocities)
    vel_i1 = velocity(1,2:numXCells+2);
    
    %create pressure vector with ghost cells
    pressure = [pressure(1,1) pressure pressure(1,size(pressure,2))];
    %copy the pressure in every i-th cell 
    %(i.e., the left cell pressures)
    press_i = pressure(1,1:numXCells+1);
    %copy the pressure in every i+1-th cell 
    %(i.e., the right cell pressures)
    press_i1 = pressure(1,2:numXCells+2);
    
    %use density, velocity, and pressure vectors 
    %to determine the energy in every cell
    energy = 0.5*currDens.*velocity.^2+pressure./(par.gamma-1);
    %copy the energies in every i-th cell 
    %(i.e., the left cell energies)
    en_i = energy(1,1:numXCells+1);
    %copy the energies in every i+1-th cell 
    %(i.e., the right cell energies)
    en_i1 = energy(1,2:numXCells+2);
    
    
    %use the density, pressure, and energy vectors 
    %to determine the enthalpy in every cell
    enthalpy = (energy+pressure)./currDens;
    %copy the enthalpies in every i-th cell 
    %(i.e., the left cell enthalpies)
    entha_i = enthalpy(1,1:numXCells+1);
    %copy the enthalpies in every i+1-th cell 
    %(i.e., the right cell enthalpies)
    entha_i1 = enthalpy(1,2:numXCells+2);
    
    %compute the roe-averaged velocity across every cell boundary
    roe_vel = (sqrt(dens_i).*vel_i+sqrt(dens_i1).*vel_i1)./(sqrt(dens_i)+sqrt(dens_i1));
    %compute the roe-averaged enthalpy across every cell boundary
    roe_entha = (sqrt(dens_i).*entha_i+sqrt(dens_i1).*entha_i1)./(sqrt(dens_i)+sqrt(dens_i1));
    %compute the roe-averaged sound speed across every cell boundary
    roe_sound = sqrt((par.gamma-1)*(roe_entha-0.5*(roe_vel).^2));
    
    %compute the eigenvalues for every inter-cell interface
    lambda_1 = abs(roe_vel-roe_sound);
    lambda_2 = abs(roe_vel);
    lambda_3 = abs(roe_vel+roe_sound);
    %apply entropy fix (with an epsilon of 0.5)
    epsilon =0.5;
    %fix lambda 1
    lambda_epsilon1 = lambda_1-epsilon;
    index_1 = lambda_epsilon1<0;
    lambda_1(index_1) = epsilon+lambda_1(index_1).^2/(2*epsilon);
    %fix lambda 2
    lambda_epsilon2 = lambda_2-epsilon;
    index_2 = lambda_epsilon2<0;
    lambda_2(index_2) = epsilon+lambda_2(index_2).^2/(2*epsilon);
    %fix lambda 3
    lambda_epsilon3 = lambda_3-epsilon;
    index_3 = lambda_epsilon3<0;
    lambda_3(index_3) = epsilon+lambda_3(index_3).^2/(2*epsilon);
    
    %create the eigenvectors for every inter-cell interface
    del_U1 = dens_i1-dens_i;
    del_U2 = dens_i1.*vel_i1-dens_i.*vel_i;
    del_U3 = en_i1 - en_i;
    
    %compute the wave amplitudes
    alpha2 = (par.gamma-1)./((roe_sound).^2).*(del_U1.*(roe_entha-roe_vel.^2)+roe_vel.*del_U2-del_U3);
    alpha1 = 1./(2*roe_sound).*(del_U1.*(roe_vel+roe_sound)-del_U2-roe_sound.*alpha2);
    alpha3 = del_U1 - (alpha1+alpha2);
    %take the sum of the eigenvectors multiplied by their 
    %respective eigenvalues and amplitudes
    K_1 = [ones([1,size(dens_i,2)]);roe_vel-roe_sound;roe_entha-roe_vel.*roe_sound];
    K_2 = [ones([1,size(dens_i,2)]);roe_vel;0.5*(roe_vel).^2];
    K_3 = [ones([1,size(dens_i,2)]);roe_vel+roe_sound;roe_entha+roe_vel.*roe_sound];
    Sum = 0.5*(alpha1.*K_1.*lambda_1+alpha2.*K_2.*lambda_2+alpha3.*K_3.*lambda_3);
    %create the flux vector F_L for every left cell
    %(should be a 3x(numXCells+1) matrix)
    F_L = [dens_i.*vel_i;dens_i.*(vel_i).^2+press_i;vel_i.*(en_i+press_i)];
    
    %create the flux vector F_R for every right cell
    %(should be a 3x(numXCells+1) matrix)
    F_R = [dens_i1.*vel_i1;dens_i1.*(vel_i1).^2+press_i1;vel_i1.*(en_i1+press_i1)];
    
    %calculate the flux through every inter-cell interface
    F_in = 0.5*(F_L+F_R)-Sum;
    
    %get the inter-cell fluxes for the left boundaries    
    F_LB = F_in(:,1:end-1);
    %get the inter-cell fluxes for the right boundaries 
    F_RB = F_in(:,2:end);
    
    %determine the max wavespeed
    max1 = max(lambda_1);
    max2 = max(lambda_2);
    max3 = max(lambda_3);
    newmax = [max1,max2,max3];
    lambdamax = max(newmax);
    %calculate the timestep size
    
    dt = par.cfl*dx/lambdamax;
    
    %use Godunov's formula to update the conservative variables
    U = U+dt./dx.*(F_LB-F_RB);

    %derive the new primitive variable values from 
    %the updated conservative variable matrix
    currDens = U(1,:);
    velocity = U(2,:)./currDens;
    pressure = (U(3,:)-0.5*currDens.*velocity.^2).*(par.gamma-1);
    
    %update the time
    t =t+dt;
    
    %update the cycle number
    cycle = cycle +1;
      
end

%After the simulation terminates, output the final 
%density, velocity, and pressure vectors
density=currDens;
velocity = velocity;
pressure = pressure;
disp(cycle)
disp(t)
end

