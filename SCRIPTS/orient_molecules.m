function sm = orient_molecules(sm, n_mol)
% NAME:
%	orient_molecules
%
% PURPOSE:
%	Orient SMs
%
% INPUTS:
%	sm: the SMs
%   n_mol: # of SMs
%
% OUTPUTS:
%	sm = the updated SMs
%
% MODIFICATION HISTORY:
%	D.Bourgeois, May 2012.

%distribute theta angles in degrees according to polar coordinates (range [-90° 90°])
%see for exemple: http://en.wikibooks.org/wiki/Mathematica/Uniform_Spherical_Distribution
theta= 180/pi*2*asin(sqrt(rand(1,n_mol)))-90; % theta is with respect to sample plane

%uniformly distribute phi angles in degrees (range [-180° 180°])
phi=180*(2*rand(1,n_mol)-1); 

for i=1:n_mol
    sm(i).theta=theta(i);
    sm(i).phi=phi(i);
end

