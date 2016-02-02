%------------------------------
% GEOL 5700 
% Model the growth/drowning of coral platforms on a subsiding volcano
%------------------------------
% Franklin Hinckley
%------------------------------

%% Clean up workspace
clearvars
close all
clc

%% Define volcanic cone
posVec = -1000:1:1000;
volcDepth = -1e-7*posVec.^4 + 3e-3*posVec.^2 + 50;
plot(posVec(500:1500),volcDepth(500:1500))
ylim([-100 100])
subRate = 1/1000; % subsidence rate [m/yr]

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
tSim = 100000; % simulation duration [yr]
tVec = 0:tStep:tSim;

%% Set up initial conditions
carbThick = zeros(size(posVec));

%% Run simulation
% Initialize output matrices
carbThick_save = zeros(length(posVec),length(tVec));
SL_save = zeros(length(tVec),1);
volcDepth_save = zeros(length(posVec),length(tVec));

% Run simulation
for ii = 1:length(tVec)            
    % Evaluate sea level 
    SL = SL_func(tVec(ii));
    
    % Evaluate coral growth rate
    G = zeros(size(carbThick));
    for jj = 1:length(G)
        zCarb = -volcDepth(jj) + SL - carbThick(jj);
        if zCarb <= 0 % check if coral is above the water
            G(jj) = 0;
        else
            G(jj) = coralGrowth(zCarb);
        end
    end
    
    % Determine new platform thickness
    carbThick = carbThick + G*tStep;
                
    % Subsidence of volcano
    volcDepth = volcDepth - subRate*tStep;
    
    % Save results
    SL_save(ii) = SL;
    carbThick_save(:,ii) = carbThick;
    volcDepth_save(:,ii) = volcDepth;
    
end

%% Plots
% Animation
pL = 500; % low plot index
pH = 1500; % high plot index
figure
for ii = 1:length(tVec)
    plot(posVec(pL:pH),carbThick_save(pL:pH,ii)+volcDepth_save(pL:pH,ii),'r')
    hold on
    plot(posVec(pL:pH),SL_save(ii)*ones(size(posVec(pL:pH))),'-b')
    plot(posVec(pL:pH),volcDepth_save(pL:pH,ii),'-k')
    hold off
    xlim([posVec(pL) posVec(pH)])
    ylim([-500 A+10])
    title('Carbonate Platforms','Fontsize',14)
    ylabel('Elevation [m]','Fontsize',12)
    xlabel('Distance from Center [m]','Fontsize',12)
%    M(:,ii) = getframe(gcf);
    pause(0.01)    
end