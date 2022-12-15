function sm_par=update_photophysics_pct(sm, sm_par, im_par)

%
% PURPOSE:
%   Update the photophysical status for the different dyes
%
% INPUTS:
%   sm: the single molecules
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   sm_par: the single molecules parameters updated for photophysical Status
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.
%	D.Bourgeois, September 2022, work with cell array.

bleached_idx=39;
activated_idx=40;
blinked_idx=41;

frame=im_par.current_frame;
sm_par.photophysical_status.n_bleached(frame)=sum([sm{bleached_idx,:}]==frame); %  # of bleached dye molecules in frame
sm_par.photophysical_status.n_blinked(frame)=sum([sm{blinked_idx,:}]==1); %  # of bleached dye molecules in frame
sm_par.photophysical_status.n_activated(frame)=sum([sm{activated_idx,:}]==frame); %  # of bleached dye molecules in frame
sm_par.photophysical_status.cum_n_bleached(frame)=sum([sm{bleached_idx,:}]~=0); %  cumulated # of bleached dye molecules in frame
sm_par.photophysical_status.cum_n_activated(frame)=sum([sm{activated_idx,:}]~=0); %  cumulated # of activated dye molecules in frame
disp(['# of bleached molecules in this frame (',sm_par.fluorophore_name, '): ', num2str(sm_par.photophysical_status.n_bleached(frame))]);
disp(['# of activated molecules in this frame (',sm_par.fluorophore_name, '): ', num2str(sm_par.photophysical_status.n_activated(frame))]);
disp(['# of blinked molecules in this frame (',sm_par.fluorophore_name, '): ', num2str(sm_par.photophysical_status.n_blinked(frame))]);
disp(['Cumulated # of bleached molecules (',sm_par.fluorophore_name, '): ', num2str(sm_par.photophysical_status.cum_n_bleached(frame))]);
disp(['Cumulated # of activated molecules (',sm_par.fluorophore_name, '): ', num2str(sm_par.photophysical_status.cum_n_activated(frame))]);
