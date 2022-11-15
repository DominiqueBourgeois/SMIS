function [acc_dose, lasers]=update_lasers_accumulated_dose(lasers, im_par)

% PURPOSE:
%   Get accumulated dose delivered by lasers
%
% INPUTS:
%	im_par: the imaging parameters
%   lasers: the lasers
%
% OUTPUTS:
%   lasers: the updated lasers for accumulated dose at the end of the current frame
%   acc_dose: the spatially dependent accumulated dose at the beginning of the current frame
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.
%	D.Bourgeois, June 2020. Following bug correction in get_beam_profile,
%	update calculation of frametime_dose and addtime_dose.


n_lasers=size(lasers,2);
% first evaluate the dose delivered during frame (frametime + addtime)
frametime_dose=zeros(1,n_lasers);
addtime_dose=zeros(1,n_lasers);
for k=1:n_lasers % add contributions from all lasers
    % exposure time to the laser in current frame [s]
    frametime_exposure=1e-3*(lasers(k).on_during_frametime*lasers(k).frametime_duration);
    addtime_exposure=1e-3*(lasers(k).on_during_addtime*lasers(k).addtime_duration);
    
    % dose [J] at peak maximum of beam in frame
    frametime_dose(k)=frametime_dose(k)+frametime_exposure*lasers(k).max_beam_profile_watt*...
        lasers(k).sequence(im_par.current_frame)/100;
    addtime_dose(k)=addtime_dose(k)+addtime_exposure*lasers(k).max_beam_profile_watt*...
        lasers(k).sequence(im_par.current_frame)/100;
end

% the total accumulated dose delivered for all lasers at this point of the stack acquisition is
acc_dose=zeros(im_par.n,im_par.m); % this is spatially dependent
for k=1:n_lasers % add contributions from all lasers
    acc_dose=acc_dose+lasers(k).accumulated_dose.*lasers(k).beam_profile/lasers(k).max_beam_profile;
end
% finally update the accumulated dose
for k=1:n_lasers % add contributions from all lasers
    lasers(k).accumulated_dose=lasers(k).accumulated_dose+frametime_dose(k)+addtime_dose(k);
end
