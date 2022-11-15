function [sms, sm_par] = process_diffusion_pct(sms, sm_par, im_par, display_par)

%
% PURPOSE:
%   Main script to process diffusion state from single molecules
%
% INPUTS:
%   sms: the single molecules 
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%   display_par: the smis display parameters
%
% OUTPUTS:
%   sms: the single molecules updated for diffusion state
%   sm_par: the sm parameters updated for field w_act
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2022, optimized for parallel computing

% Open figure if necessary
if display_par.show_diff_image==1
    figure(display_par.diff_figure.Number);
end

n_fluorophores=numel(sms);

%Update status of activated molecules 
for i=1:n_fluorophores
    %Extract the useful indices
    w_idx = find(matches(sm_par(i).sm_fn,{'bleached','activated'})==1);
    bleached_idx=w_idx(1);
    activated_idx=w_idx(2);

    if im_par.move_non_activated_molecules==1
        sm_par(i).w_act = find([sms(i).sm_cell{bleached_idx,:}]==0 | [sms(i).sm_cell{bleached_idx,:}]==im_par.current_frame); % Also treat those fluorophores that just bleached
    else
        sm_par(i).w_act = find([sms(i).sm_cell{activated_idx,:}] & ([sms(i).sm_cell{bleached_idx,:}]==0 | [sms(i).sm_cell{bleached_idx,:}]==im_par.current_frame)); % Also treat those fluorophores that just bleached
    end
end


% Process the diffusing molecules
if im_par.simul_3D==0 % Set to 1 if diffusion is confined within pattern
    sms = move_diffusing_sm_2D_pct(n_fluorophores, sms, sm_par, im_par, display_par);
else
    sms = move_diffusing_sm_3D_pct(n_fluorophores, sms, sm_par, im_par, display_par);
end


if im_par.use_diffuse_psf==1 % Correct for xyz_track positions to the mean position instead of the end position during frame time
    sms=update_diffuse_track_pct(n_fluorophores, sms, sm_par, im_par);
end






