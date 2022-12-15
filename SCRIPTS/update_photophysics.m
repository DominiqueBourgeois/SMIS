function sm_par=update_photophysics(sms, sm_par, im_par)

%
% PURPOSE:
%   Update the photophysical status for the different dyes
%
% INPUTS:
%   sms: the single molecules (with coordinates on high-resolution image) in raster units
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   sm_par: the single molecules parameters updated for photophysical Status
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.

n_dyes=size(sm_par,2);
frame=im_par.current_frame;
for k=1:n_dyes
    sm=sms(k).sm;
    sm_par(k).photophysical_status.n_bleached(frame)=sum([sm.bleached]==frame); %  # of bleached dye molecules in frame
    sm_par(k).photophysical_status.n_blinked(frame)=sum([sm.blinked]==1); %  # of bleached dye molecules in frame
    sm_par(k).photophysical_status.n_activated(frame)=sum([sm.activated]==frame); %  # of bleached dye molecules in frame
    sm_par(k).photophysical_status.cum_n_bleached(frame)=sum([sm.bleached]~=0); %  cumulated # of bleached dye molecules in frame
    sm_par(k).photophysical_status.cum_n_activated(frame)=sum([sm.activated]~=0); %  cumulated # of activated dye molecules in frame
    disp(['# of bleached molecules in this frame (',sm_par(k).fluorophore_name, '): ', num2str(sm_par(k).photophysical_status.n_bleached(frame))]);
    disp(['# of activated molecules in this frame (',sm_par(k).fluorophore_name, '): ', num2str(sm_par(k).photophysical_status.n_activated(frame))]);
    disp(['# of blinked molecules in this frame (',sm_par(k).fluorophore_name, '): ', num2str(sm_par(k).photophysical_status.n_blinked(frame))]);
    disp(['Cumulated # of bleached molecules (',sm_par(k).fluorophore_name, '): ', num2str(sm_par(k).photophysical_status.cum_n_bleached(frame))]);
    disp(['Cumulated # of activated molecules (',sm_par(k).fluorophore_name, '): ', num2str(sm_par(k).photophysical_status.cum_n_activated(frame))]);
end