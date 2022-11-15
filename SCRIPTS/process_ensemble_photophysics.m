function [ens, sm_par] = process_ensemble_photophysics(ens, lasers, sm_par, im_par)

%
% PURPOSE:
%   Main script to process photophysical state from single molecules
%
% INPUTS:
%   ens: the ensemble data with state populations
%   lasers: the lasers
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   P: the updated state populations
%	sm_par: the sm parameters eventually updated for sampling rate
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2021.


% Each sub population will be processed individually
for k=1:sm_par.n_subpop
    if im_par.addtime>0 %Process photophysics during addtime
        im_par.during_frametime=0; 
        [ens.sp(k), sm_par]=get_ensemble_state_evolution(ens.sp(k),lasers, sm_par, im_par);
    end
    %Now process photophysics during frametime
    im_par.during_frametime=1;
    [ens.sp(k),sm_par]=get_ensemble_state_evolution(ens.sp(k), lasers, sm_par, im_par);
end




