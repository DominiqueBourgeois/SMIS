function  sm = update_tracks_bleached_sm_3D(sm, im_par)

% PURPOSE:
%	Update tracks for molecules bleached in current frame
%
% INPUTS:
%   sm: the single molecules (with coordinates on high-resolution image)
%	im_par: the imaging parameters
%
% OUTPUTS:
%   sm: the single molecules updated for: x_track, y_track
%
% MODIFICATION HISTORY:
%	D.Bourgeois, October 2020: version > simulate_palm_vsn16.2

for k=1:numel(sm) % Go for all molecules   
    sm(k).x_track=vertcat(sm(k).x_track, [sm(k).x,im_par.current_frame]);
    sm(k).y_track=vertcat(sm(k).y_track, [sm(k).y,im_par.current_frame]);
    sm(k).z_track=vertcat(sm(k).z_track, [sm(k).z,im_par.current_frame]);
end

