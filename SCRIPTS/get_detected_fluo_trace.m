function ft_det=get_detected_fluo_trace(ft,state, frametime,min_detectable_on)

% PURPOSE:
%   Get detected fluorescence trace for a certain state
%
% INPUTS:
%	ft: initial fluorescence trace
%	state: the state id
%   frametime: frametime in [s]
%   min_detectable_on: minimum fraction of time the state should be present during frametime to be considered detectable
%
% OUTPUTS:
%	ft_det: the detected trace for the state in question
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2020.


ft_det=ft;

t_start=ft(1,1); % Start of trace [s]
t_end=ft(1,end);

f_start=1+fix(t_start/frametime); % First frame
f_end=1+fix(t_end/frametime);

% Go through all frames during which state appears
for k=f_start:f_end
    t1=(k-1)*frametime;
    t2=k*frametime;
    ft_det=get_main_state(ft_det,state, t1,t2,min_detectable_on*frametime);
end


