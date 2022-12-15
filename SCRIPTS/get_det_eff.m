function [a, b, mean_eff] = get_det_eff(opening_angle, mic_transmission)
% NAME:
%   get_det_eff
%
% PURPOSE:
%   This program calculates detector efficiency based on Fourkas, Opt. Lett. 2001 
%   for unpolarized detection
%
% INPUTS:
%   opening_angle: half opening angle of objective
%   mic_transmission: efficiency of the detection system (% of photons
%       captures by the objective that reach the EMCCD detector)
%
% OUTPUTS:
%	a: efficiency of collection
%
% MODIFICATION HISTORY:
%	D.Bourgeois, May 2012.
%

a=1/6-1/4*cos(opening_angle)+1/12*(cos(opening_angle))^3;
b=1/8*cos(opening_angle)-1/8*(cos(opening_angle))^3;
mean_eff=mic_transmission*2*(a+1/3*b); % the mean of <sin(theta)^2> over spherical coordinates is 1/3
end