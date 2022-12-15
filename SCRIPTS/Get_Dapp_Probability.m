function p = Get_Dapp_Probability(Dapp,D,tracklength)

% PURPOSE:
% Get probability of measuring an apparent diffusion coeff Dapp, given a
% given D and a tracklength, when Dapp is measured from the mean jump
% distance square <JD²>
%
% INPUTS:
%   Dapp: The apparent diffusion coeff
%   D: The true diffusion coeff
%   tracklength: the tracklength [frames]
%
% OUTPUTS:
%	p: the probability of measuring Dapp
%
% MODIFICATION HISTORY:
%	D.Bourgeois, January 2022.

if tracklength>1
    p=1/factorial(tracklength-1)*(tracklength/D)^tracklength*Dapp^(tracklength-1)*exp(-tracklength*Dapp/D);
else
    p=nan;
end
end