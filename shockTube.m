%shockTube.m
%main program for Assignment 3 shock tube simulation

%shockTube.m
%main program for Assignment 3 shock tube simulation

function shockTube

%Retrieve the shock tube parameters from the setup function.
par = setup();

%Use the functions exactSolution.m and roeSolution.m to find, respectively,
%the exact and simulated results for the state of the shock tube at the
%specified time.

[exactDensity, exactVelocity, exactPressure] = exactSolution(par);

[simDensity, simVelocity, simPressure] = roeSolution(par);

%Use plots to compare the exact and simulated results.
xCoords = par.cellCoords;
a1=plot(xCoords,exactDensity);M1='exact';
%a1=plot(xCoords,exactVelocity);M1='exact';
%a1=plot(xCoords,exactPressure);M1='exact';
hold on;
a2=plot(xCoords,simDensity);M2='roe';
%a2=plot(xCoords,simVelocity);M2='roe';
%a2=plot(xCoords,simPressure);M2='roe';
%hold off;
title('Sod''s shock Tude at t=0.2(exact vs roe solution)')
xlabel('x')
ylabel('Density')
%ylabel('Pressure')
%ylabel('Velocity')
legend([a1,a2],M1,M2)
%legend('t=0.2','t=0.3','t=0.4','t=0.5','t=0.6')
end

function par = setup()
    %Set the initial shock tube state here.
    par.densL = 1.0; %left chamber density
    par.vxL  = 0.75; %left chamber velocity
    par.presL = 1.0; %left chamber pressure

    par.densR = 0.125; %right chamber density
    par.vxR  = 0.0; %right chamber velocity
    par.presR = 0.1; %right chamber pressure

    %Set gas constant.
    par.gamma = 1.4; %specific heat ratio
    
    %It is also convenient to calcuate the initial
    %sound speed in the left chamber.
    par.csL = sqrt(par.gamma*par.presL/par.densL);

    %Set simulation parameters.
    par.cfl = 0.5; %Courant-Friedrichs-Lewy number

    par.maxCycles = 1000; %maximum number of cycles/iterations

    par.numXCells = 400; %number of cells in the simulation grid
    par.xMin = 0.0; %min x-position 
    par.xMax = 1.0; %max x-position
    par.dx = (par.xMax-par.xMin)/(par.numXCells-1); %cell width
    %Create a vector that containes the x-positions of every 
    %cell in the simulation grid.
    par.cellCoords = par.xMin:par.dx:par.xMax;

    %Set the diaphragm position.
    par.x0 = 0.5;

    %Choose the time at which you want to evaluate the 
    %shock tube solution. 
    par.t = 0.2;
end

