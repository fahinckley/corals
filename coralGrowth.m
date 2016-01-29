%------------------------------
% Determine the growth rate of the coral platform at the specified depth
%------------------------------
% Franklin Hinckley
% 29 January 2016
%------------------------------
% Inputs:
%   z (n x 1 double array): depths [m]
% Outputs:
%   G (n x 1 double array): Growth rate at each depth [m/yr]
%------------------------------

function [G] = coralGrowth(z)

%% Define parameters
Gm = 7.5/1000; % maximum growth rate [m/yr]
k = 0.15; % extinction coefficient [1/m]
I0 = 2000; % surface light intensity [uE/m^2/s]
Ik = 300; % saturating light intensity [uE/m^2/s]

%% Compute growth rate
G = Gm.*tanh(I0.*exp(-k.*z)./Ik);
