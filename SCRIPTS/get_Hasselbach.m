function [A, AH]=get_Hasselbach(pH, pKa, Hill)

% PURPOSE:
%       Get concentration of anionic and neutral species as function of pKa
%       and pH
% INPUTS:
%   pH, pKa, Hill coeff
%
% OUTPUTS:
%	A, AH : fractional concentration of anionic and neutral
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2020
%-

r=10.^((pH-pKa).*Hill);
AH=1./(1+r);
A=1-AH;
