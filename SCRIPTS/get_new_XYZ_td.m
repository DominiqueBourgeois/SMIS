function sm_td=get_new_XYZ_td(sm, sm_td, im_par, pair_id)

% PURPOSE:
% Get new position in 3D for a tandem molecule linked to a fluorophore that was moved by diffusion
%
% INPUTS:
%	sm: the already moved fluorophore
%	sm_td: the tandem fluorophore that needs to be moved
%   im_par: imaging parameter
%   pair_id: the index of fluorophore_pairs corresponding to the reference and target SMs (index in sm_par)
%
% OUTPUTS:
%	sm_td: the updated tandem fluorophore that needs to be moved
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15


r1=im_par.d1d2_dist(pair_id)*im_par.binning/im_par.raster; %in pixels for the unbinned image
r2=im_par.d1d2_dist_sig(pair_id)*im_par.binning/im_par.raster; %in pixels for the unbinned image
d = r1 + r2*randn(1,1); % define random distances
theta = 2*pi*rand(1,1); % define random angle
phi = -pi/2+pi*rand(1,1); % define random angle
sm_td.x = sm.x + d*cos(theta)*cos(phi);
sm_td.y = sm.y + d*sin(theta)*cos(phi);
sm_td.z = sm.z + d*sin(phi);
% do the same for intermediate xy if necessary
if im_par.use_diffuse_psf==1
    sm_td.sub_x = sm.sub_x + d*cos(theta)*cos(phi);
    sm_td.sub_y = sm.sub_y + d*sin(theta)*cos(phi);
    sm_td.sub_z = sm.sub_z + d*sin(phi);
end
end
