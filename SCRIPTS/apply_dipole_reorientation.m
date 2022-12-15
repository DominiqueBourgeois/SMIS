function [sms] = apply_dipole_reorientation(sms, sm_par, im_par)

% PURPOSE:
%	reorient dipole molecules according to their probability to reorient
%
% INPUTS:
%   sms: the single molecules (with current theta)
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   sms: the single molecules updated for: theta
%
% MODIFICATION HISTORY:
%	Jip Wulffele, 12/17/2020
%	Dominique Bourgeois, 03/11/2021
%
%Strategy:
%=> for each molecule decide wheter or not to change orientation
%=> if reorient=yes => chance sm.theta 

%only deal with activated molecules
w_act=sm_par.w_act;

k = sm_par.jump_rate; % [-1] Rate of reorientation
T=(im_par.frametime + im_par.addtime)/1000; % Total frame time [s]
sampling_rate = max([10*k,1/T]); % Sampling rate: Try at least once per frame, and 10 times faster than the jump rate

p = k/sampling_rate; % Probability of reorientation for each trial

N_trials=round(T*sampling_rate); % Number of trials during frame time: we will do a test every 1/sampling_rate [s].

%loop over active molecules
for j=1:length(w_act)  
    Test=rand(1,N_trials); % series of random numbers
    r = find(Test <= p); % Results of the test
    
    if ~isempty(r) % Reorientation occurred
        currenttheta=sms.sm(w_act(j)).theta;
        
        newthetas= 180/pi*2*asin(sqrt(rand(1,numel(r))))-90; % theta is with respect to sample plane
        
        fullthetas= [currenttheta,newthetas]; % The new angles distribute between -90 and 90°
        t_jump=sort((r))*T/N_trials; % The reorientation times
        delta_t=circshift(t_jump,-1,2)-t_jump;
        delta_t=[t_jump(1),delta_t(1:end-1),T-t_jump(end)];
        sms.sm(w_act(j)).theta=sum(abs(fullthetas)*delta_t)/T; % Weighted-average of the absolute theta values (mean of theta and -theta should be theta, not zeo !      
    end
end
    
  