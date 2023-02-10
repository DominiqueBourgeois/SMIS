function sm_par = check_sampling_rate(sms, lasers, sm_par, im_par)

%
% PURPOSE:
%   Check sampling rate of sm
%
% INPUTS:
%   sms: the single molecules (with coordinates on high-resolution image) in raster units
%   lasers: the lasers
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%	sm_par: the sm parameters updated for sampling rate
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2022.


% Define an average molecule
sm0=sms.sm(1); 
sm0.x=mean([sms.sm.x]); 
sm0.y=mean([sms.sm.y]); 
if im_par.simul_3D==1
    sm0.z=mean([sms.sm.z]);
end
if sm_par.anisotropy==1 % take anisotropy into account
    sm0.theta=mean([sms.sm.theta]);
    sm0.phi=mean([sms.sm.phi]);
end

im_par.current_frame=1; % Get current frametime to 1 for the get_number_of_absorbed_photons.m script
% Get the sampling rates
if im_par.addtime>0 % Get initial sampling rate for addime
    im_par.during_frametime=0;
    sm_par=get_sampling_rate(sm0,lasers, sm_par, im_par);
end
% Get initial sampling rate for frametime
im_par.during_frametime=1;
sm_par=get_sampling_rate(sm0, lasers, sm_par, im_par);



