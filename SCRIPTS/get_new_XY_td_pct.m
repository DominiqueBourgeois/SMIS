function sm_td=get_new_XY_td_pct(sm, sm_td, im_par, pair_id)

% PURPOSE:
% Get new position for a tandem molecule linked to a fluorophore that was moved by diffusion
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
%	D.Bourgeois, September 2022, optimized for parallel computing

x_idx=1;
y_idx=2;
sub_x_idx=3;
sub_y_idx=4;

r1=im_par.d1d2_dist(pair_id)*im_par.binning/im_par.raster; %in pixels for the unbinned image
r2=im_par.d1d2_dist_sig(pair_id)*im_par.binning/im_par.raster; %in pixels for the unbinned image
d = r1 + r2*randn(1,1); % define random distances
theta = 2*pi*rand(1,1); % define random angle
sm_td{x_idx} = sm{x_idx} + d.*cos(theta);
sm_td{y_idx} = sm{y_idx} + d.*cos(theta);
% do the same for intermediate xy if necessary
if im_par.use_diffuse_psf==1
    sm_td{sub_x_idx} = sm{sub_x_idx} + d.*cos(theta);
    sm_td{sub_y_idx} = sm{sub_y_idx} + d.*cos(theta);
end
