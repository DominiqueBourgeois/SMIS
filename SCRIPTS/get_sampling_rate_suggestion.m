function [suggested_sampling_rate,suggested_frametime] = get_sampling_rate_suggestion(frametime,min_sampling_points)

% PURPOSE:
% Get adapted sampling rate with the following constraints:
% Sampling rate must be an integer
% Number of samples during frame time must be at least 100
%
% INPUTS:
% frametime : The frametime or addtime [ms]
% min_sampling_points: minimum number of sampling points per frametime or addtime
%
% OUTPUTS:
%	The suggested sampling rate [Hz]
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2021

suggested_frametime=[];
k=min_sampling_points;
while 1e3*k/frametime~=floor(1e3*k/frametime) && 1e3*k/frametime<1e+7 % Do not accept a sampling rate more than 10 MHz
    k=k+1;
end
suggested_sampling_rate=1e3*k/frametime;

% In that case, the frame time should be adjusted but we do not want a too small adjustment 
% (only two digits are allowed and the second digit should be an even digit)
if suggested_sampling_rate>=1e7    
    if rem(floor(100*frametime),2)==0
        suggested_frametime=0.01*floor(100*frametime);
    else
        suggested_frametime=0.01*floor(100*frametime+1);
    end
end



