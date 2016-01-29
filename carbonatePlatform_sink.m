%------------------------------
% GEOL 5700 
% Model the growth/drowning of coral platforms
%------------------------------
% Franklin Hinckley
%------------------------------

%% Clean up workspace
clearvars
close all
clc

%% Define plate motion
w0 = 10*1000; % deflection at load point [m]
conRate = 8/100; % convergence rate of load [m/yr]
loadPos = 100*1000; % distance from load to land edge [m]

%% Define function for sea level
A = 120/2; % half-amplitude of sea level variation [m]
P = 20000; % period of oscillation [yr]
B = (2*pi)/20000; 
C = 0; % phase of oscillation at initial time [rad]
D = 0; % mean sea level relative to reference height [m]
SL_func = @(t) A * sin(B * t + C) + D;

%% Define simulation parameters
% Simulation time
tStep = 100; % simulation time step [yr] 
tSim = 500000; % simulation duration [yr]
tVec = 0:tStep:tSim;

%% Set up initial conditions
posVec = fliplr(0:conRate*tStep:loadPos);
carbThick = zeros(size(posVec));

%% Run simulation
% Initialize output matrices
carbThick_save = zeros(length(posVec),length(tVec));
SL_save = zeros(length(tVec),1);

% Evaluate plate depth
plateDepth = plateFlexure(w0,posVec);

% Run simulation
for ii = 1:length(tVec)            
    % Evaluate sea level 
    SL = SL_func(tVec(ii));
    
    % Evaluate coral growth rate
    G = zeros(size(carbThick));
    for jj = 1:length(G)
        zCarb = plateDepth(jj) + SL - carbThick(jj);
        if zCarb <= 0 % check if coral is above the water
            G(jj) = 0;
        else
            G(jj) = coralGrowth(zCarb);
        end
    end
    
    % Determine new platform thickness
    carbThickNew = carbThick + G*tStep;
    
    % Remove points under the load and add new points to the land edge
    carbThick = zeros(size(carbThickNew));
    carbThick(posVec > 0) = carbThickNew(posVec > 0);
            
    % Save results
    SL_save(ii) = SL;
    carbThick_save(:,ii) = carbThick;
    
    % Shift points for movement towards load point
    carbThick = [0 carbThick(1:end-1)];
    
end

%% Plots
% Final state
figure
plot(posVec,carbThick_save(:,end))

% Animation
pL = 1500; % low plot index
pH = 3000; % high plot index
figure
for ii = 1:length(tVec)
    plot(posVec(pL:pH)/1000,carbThick_save(pL:pH,ii)-plateDepth(pL:pH)','r')
    hold on
    plot(posVec(pL:pH)/1000,SL_save(ii)*ones(size(posVec(pL:pH))),'-b')
    plot(posVec(pL:pH)/1000,-plateDepth(pL:pH),'-k')
    hold off
    xlim([posVec(pH)/1000 posVec(pL)/1000])
    ylim([-500 A+10])
    title('Carbonate Platforms','Fontsize',14)
    ylabel('Elevation [m]','Fontsize',12)
    xlabel('Distance from Load [km]','Fontsize',12)
%    M(:,ii) = getframe(gcf);
    pause(0.01)    
end