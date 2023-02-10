function [sm_cell,det_im,ok]=get_frame_image_pct(sm_cell, sm_pattern_indices, sm_par, det_im, im_par)

%
% PURPOSE:
%   Get images as recorded by the detector
%
% INPUTS:
%   sm: the single molecules 
%   sm_pattern_indices: indices of virtual sample subpatterns
%	sm_par: the sm parameters
%	det_im: the detector images
%	im_par: the imaging parameters
%
% OUTPUTS:
%	det_im: the updated detector images for emccd_im_ch1 (image on channel 1) and emccd_im_ch2
%   (image on channel 2)
%   sms: The updated sms
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.
%	D.Bourgeois, May 2021. Introduce sm.n_phot_det_ch1 and sm.n_phot_det_ch2 in update_image and update_sub_image;
%	D.Bourgeois, September 2022, use cell arrays instead of structures
%	D.Bourgeois, February 2023, introduce sm_pattern_indices, now disconnected from sm_par

% Note: parallel computing cannot be used here because a single image is
% modified

ok=1; % In case of error message (necessary for get_linkage_error)

%Extract the needed fields for sm_cell
sm=sm_cell(sm_par.idx.get_frame_image.indices,:);

% Get the proper indices for donor and acceptor
x_idx=1;
y_idx=2;
z_idx=3;
sub_x_idx=4;
sub_y_idx=5;
sub_z_idx=6;
n_phot_ch1_idx=9;
n_phot_ch2_idx=10;
bleached_idx=17;
% n_phot_ch1_idx=8;
% n_phot_ch2_idx=9;
% bleached_idx=16;


for k=1:sm_par.n_mol_eff
    if sm{bleached_idx,k}==0 || sm{bleached_idx,k}==im_par.current_frame
        if im_par.add_diffusion==0 || im_par.use_diffuse_psf==0 || isempty(sm{sub_x_idx,k}) % general case
            [x,y,z]=get_coordinates_on_detector_pct(sm{x_idx,k},sm{y_idx,k},sm{z_idx,k},im_par.binning);
            if any(sm{n_phot_ch1_idx,k})
                [sm(:,k),det_im,ok]=update_image_pct(x,y,z,sm(:,k), det_im, im_par, sm_pattern_indices, sm_par, 1);
                if ok==0
                    break
                end
            end
            if im_par.two_channel==1 && any(sm{n_phot_ch2_idx,k})
                [sm(:,k),det_im,ok]=update_image_pct(x,y,z,sm(:,k), det_im, im_par, sm_pattern_indices, sm_par, 2);
                if ok==0
                    break
                end
            end
        else
            [x,y,z]=get_coordinates_on_detector_pct(sm{sub_x_idx,k},sm{sub_y_idx,k},sm{sub_z_idx,k}, im_par.binning);
            if any(sm{n_phot_ch1_idx,k})
                [sm(:,k),det_im]=update_sub_image_pct(x,y,z,sm(:,k), det_im, im_par, sm_pattern_indices, sm_par, 1);
            end
            if im_par.two_channel==1 && any(sm{n_phot_ch2_idx,k})
                [sm(:,k),det_im]=update_sub_image_pct(x,y,z,sm(:,k), det_im, im_par, sm_pattern_indices, sm_par, 2);
            end
        end
    end
end

%Fill up sm_cell
sm_cell(sm_par.idx.get_frame_image.indices,:)=sm;

