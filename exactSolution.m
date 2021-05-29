%exactSolution.m
%Purpose: finds the exact solution for Sod's shock tube 
%Parameters: (passed through the par structure)
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
%            t           the time at which the function evaluates
%                        the exact solution
%Outputs: density  a vector containing the exact solution for the gas
%                  density in every cell (same size as par.cellCoords param)
%         velocity a vector containing the exact solution for the gas
%                  velocity in every cell (same size as par.cellCoords param) 
%         pressure a vector containing the exact solution for the gas
%                  pressure in every cell (same size as par.cellCoords param)

function [density,velocity,pressure] = exactSolution(par)

%create empty vectors for the density, velocity,
%and pressure outputs

density = zeros([1 size(par.cellCoords,2)]);
velocity = zeros([1 size(par.cellCoords,2)]);
pressure = zeros([1 size(par.cellCoords,2)]);

%use the initial conditions to determine the density, 
%velocity, and pressure in each of the five regions
%described in the assignment document

%Region 1
%density, velocity, and pressure equal to 
%the initial left chamber density, velocity, and pressure. 
dens_1 = par.densL;
press_1 = par.presL;
u_1 = par.vxL;
%also calculate the sound speed in Region 1
cs_1 = par.csL;
%Skip Region 2 for now (depends on x-coordinates)


%Regions 3 and 4
%use the Newton's method solver to find p-star
%(with a tolerance of 0.01 and an initial p-star guess of (p_L+p_R)/2)
%(Note: the value of p-star has to be somewhere between p_L and p_R,
%       so the average (p_L + p_R)/2 is a good place to start 
%       looking for a solution.)
pStar = 0.5*(par.presL+par.presR);
pStar = pStarSolver(pStar,0.01,par);
%use p-star to find the Region 3 density and the Region 4 density
press_3 = pStar;
press_4 = press_3;
dens_3 = dens_1 .* (press_3/press_1).^(1/par.gamma);
A = 2/(par.densR*(par.gamma+1));
B = par.presR*(par.gamma-1)/(par.gamma+1);
dens_4 = par.densR*(par.presR*(par.gamma-1)+press_4*(par.gamma+1))/(press_4*(par.gamma-1)+par.presR*(par.gamma+1));

%use p-star to find vx-star, the velocity in Regions 3 and 4
uStar = par.vxR + (pStar - par.presR)*sqrt(A/(pStar+B));
u_3 = uStar;
u_4 = u_3;

%Region 5
%Set density, velocity, and pressure equal to the
%initial right chamber density, velocity, and pressure. 
dens_5 = par.densR;
u_5 = par.vxR;
press_5 = par.presR;

%Next we need to find the start and end of each region.
%Note: Before adding data to the density, velocity, and pressure vectors,
%      use assert statements to ensure that the start and end of the
%      region you're working on both fall within the domain of the 
%      shock tube system. Including these assert statements will help
%      you answer Question 1 of this assignment. 

%ensure that the shock tube grid has at least 2 cells
assert(size(par.cellCoords, 2) > 1);

%calculate the width of each cell (could also just use par.dx parameter)
dx = par.cellCoords(1,2) - par.cellCoords(1,1);

%the size of Region 1 depends on the speed of 
%the head of the rarefaction wave.

%Copy the Region 1 density, velocity, and pressure
%into the appropriate cells in the output vectors.
x_1start = par.xMin;
vHead = u_1 - cs_1;
x_1end = par.x0 + par.t*vHead;
x_1ceil = ceil(x_1end*par.numXCells);
x1range = 1:x_1ceil;
density(x1range) = dens_1;
velocity(x1range) = u_1;
pressure(x1range) = press_1;

%the size of Region 2 depends on the speed of 
%the tail of the rarefaction wave.
cs_3 = sqrt(par.gamma*press_3/dens_3);
vTail = u_3 - cs_3;
x_2end = par.x0 + par.t*vTail;
x_2ceil = x_1ceil +1;
x_2ceil2 = ceil(x_2end*par.numXCells);
x2range = x_2ceil:x_2ceil2;



%Use the x-positions provided in par.cellCoords to find the density, velocity,
%and pressure everywhere inside Region 2. (Note: Region 2 is the only
%region where the density, velocity, and pressure are not uniform.)
density(x2range) = dens_1.*(2/(par.gamma+1)+(par.gamma-1)/...
(cs_1*(par.gamma+1)).*(u_1-(par.cellCoords(x2range)-par.x0)/par.t)).^(2/(par.gamma-1));
velocity(x2range) = (2/(par.gamma+1)).*(cs_1+0.5*(par.gamma-1)*u_1+(par.cellCoords(x2range)-par.x0)/par.t);
pressure(x2range) = press_1.*(2/(par.gamma+1)+(par.gamma-1)/(cs_1*(par.gamma+1)).*(u_1-(par.cellCoords(x2range)-par.x0)/par.t)).^(2*par.gamma/(par.gamma-1));
%the size of Region 3 depends on the speed of the right-travelling
%contact discontinuity.

%Copy the Region 3 density, velocity, and pressure
%into the output vectors.
x_3ceil = x_2ceil2+1;
Vcontact = uStar;
x_3end = par.x0+par.t*Vcontact;
x_3ceil2 = ceil(x_3end*par.numXCells);
x3range = x_3ceil:x_3ceil2;
density(x3range) = dens_3;
velocity(x3range) = u_3;
pressure(x3range) = press_3;

%the size of Region 4 depends on the speed of the shock wave.
x_4ceil = x_3ceil2+1;
cs_5 = sqrt(par.gamma*press_5/dens_5);
vShock = u_5 +cs_5*sqrt((par.gamma+1)*press_4/(2*par.gamma*press_5)+(par.gamma-1)/(2*par.gamma));
x_4end = par.x0+par.t*vShock;
x_4ceil2 = ceil(x_4end*par.numXCells);
x4range = x_4ceil:x_4ceil2;
%into the output vectors.
density(x4range) = dens_4;
velocity(x4range) = u_4;
pressure(x4range) = press_4;

%Region 5 occupies the rest of the shock tube grid
x_5ceil = x_4ceil2+1;
x_5ceil2 = par.numXCells;
x5range =x_5ceil:x_5ceil2;
density(x5range) = dens_5;
velocity(x5range) = u_5;
pressure(x5range) = press_5;
%plot(par.cellCoords,density)
%plot(par.cellCoords,velocity)
%plot(par.cellCoords,pressure)
%Copy the Region 5 density, velocity, and pressure
%into the remaining cells.


end

