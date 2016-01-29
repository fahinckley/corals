%------------------------------
% Determine the flexure of a semi-infinite elastic plate under a line load
%------------------------------
% Franklin Hinckley
% 29 January 2016
%------------------------------
% Inputs:
%   w0 (double): maximum deflection [m]
%   x (n x 1 double array): distances from the load to 
%       points of interest [m]
% Outputs:
%   w (n x 1 double array): deflections at each x point [m]
%------------------------------

function [w] = plateFlexure(w0,x)

%% Define parameters and compute flexural parameter
rhoM = 3300; % density of mantle [kg/m^3]
rhoW = 1030; % density of sea water [kg/m^3]
g = 9.81; % acceleration due to gravity [m/s^2]
D = 5e22; % flexural rigidity [Nm]
alpha = (4*D/((rhoM - rhoW)*g))^(1/4);

%% Determine deflection at the specified distances from the load
w = w0.*exp(-x/alpha).*cos(x/alpha);
