function [sm, det_im, ok] = update_image_pct(x,y,z,sm, det_im, im_par, sm_par, channel)
%
% NAME:
%	update_image
% PURPOSE:
%
% INPUTS:
%   x, y, z: coordinate of the single molecule on detector image [raster
%   units]
%   sm: the single molecule (with coordinates on high-resolution image) in raster units
%	det_im: the detector images 
%	im_par: the imaging parameters
%	sm_par: the sm parameters
%   channel: set to channel number (1 for channel 1, etc)
%
% OUTPUTS:
%	det_im: the updated detector images 
%   sm: The updated single molecule
%
% MODIFICATION HISTORY:
%	D.Bourgeois, March 2012. Modified September 2019
%	D.Bourgeois, September 2020. Only update image if photons are emitted
%	D.Bourgeois, May 2021. Introduce sm.n_phot_det_ch1 and sm.n_phot_det_ch2;
%	D.Bourgeois, July 2021. Keep track of sm.n_phot_det_ch1 and sm.n_phot_det_ch2 over whole data set;
%	D.Bourgeois, September 2022, use cell arrays instead of structures. Introduce det_im. Introduce linkage error


ok=1; % In case of error message (necessary for get_linkage_error)

% Get the proper indices 
c_sp_idx=7;
theta_idx=8;
n_phot_ch1_idx=9;
n_phot_ch2_idx=10;
n_phot_det_ch1_idx=11;
n_phot_det_ch2_idx=12;
tot_phot_det_ch1_idx=13;
tot_phot_det_ch2_idx=14;
frames_on_ch1_idx=15;
frames_on_ch2_idx=16;
lx_idx=19;
ly_idx=20;
lz_idx=21;
le_set_idx=22;

% theta_idx=7;
% n_phot_ch1_idx=8;
% n_phot_ch2_idx=9;
% n_phot_det_ch1_idx=10;
% n_phot_det_ch2_idx=11;
% tot_phot_det_ch1_idx=12;
% tot_phot_det_ch2_idx=13;
% frames_on_ch1_idx=14;
% frames_on_ch2_idx=15;
% lx_idx=18;
% ly_idx=19;
% lz_idx=20;
% le_set_idx=21;

% Set the number of detected photons
n_phot_det_ch1=0;
n_phot_det_ch2=0;

% possibly apply distortion if working in channel 2
if channel==2
    [x, y]=apply_distortion(x,y,im_par);
end

% go through all fluorescence states
for i=1:sm_par.n_fluorescent_states
    % treat channel 1 first
    if channel==1 && sm{n_phot_ch1_idx}(i)>0
        %get the number of photons actually reaching the detector
        n_phot_det_ch1=n_phot_det_ch1+get_emitted_photons(sm{n_phot_ch1_idx}(i), sm{theta_idx}, im_par, sm_par);
    elseif channel==2 && sm{n_phot_ch2_idx}(i)>0 %same for channel 2
        %get the number of photons actually reaching the detector
        n_phot_det_ch2=n_phot_det_ch2+get_emitted_photons(sm{n_phot_ch2_idx}(i), sm{theta_idx}, im_par, sm_par);
    end
end

if channel==1 && n_phot_det_ch1>0
    %update sm{n_phot_det_ch1_idx}
    sm{n_phot_det_ch1_idx}=[sm{n_phot_det_ch1_idx},n_phot_det_ch1];
    sm{frames_on_ch1_idx}=[sm{frames_on_ch1_idx},im_par.current_frame];
    
    % Get the linkage error
    if sm_par.linkage_error==1
        [x,y,z,sm([c_sp_idx,lx_idx,ly_idx,lz_idx,le_set_idx]),ok] = ...
            get_linkage_error_pct(x,y,z,sm([c_sp_idx,lx_idx,ly_idx,lz_idx,le_set_idx]), im_par, sm_par);
    end

    %get the number of photons actually reaching the detector
    if im_par.simul_3D==0
        peak=get_2Dpeak(x,y,n_phot_det_ch1 ,sm_par,im_par,i,channel);
    else
        peak=get_3Dpeak(x,y,z, n_phot_det_ch1, sm_par,im_par,i,channel);
    end
    %Position this spot into the whole image
    det_im.emccd_im_ch1=position_peak(peak,[round(x),round(y)],det_im.emccd_im_ch1);
elseif channel==2 && n_phot_det_ch2>0 %same for channel 2
    %update sm{n_phot_det_ch2_idx}
    sm{n_phot_det_ch2_idx}=[sm{n_phot_det_ch2_idx},n_phot_det_ch2];
    sm{frames_on_ch2_idx}=[sm{frames_on_ch2_idx},im_par.current_frame];
    
    % Get the linkage error
    if sm_par.linkage_error==1
        [x,y,z,sm([c_sp_idx,lx_idx,ly_idx,lz_idx,le_set_idx]),ok] = ...
            get_linkage_error_pct(x,y,z,sm([c_sp_idx,lx_idx,ly_idx,lz_idx,le_set_idx]), im_par, sm_par);
    end

    %get the number of photons actually reaching the detector
    if im_par.simul_3D==0
        peak=get_2Dpeak(x,y,n_phot_det_ch2,sm_par,im_par,i,channel);
    else
        peak=get_3Dpeak(x,y,z,n_phot_det_ch2,sm_par,im_par,i,channel);
    end
    %Position this spot into the whole image
    det_im.emccd_im_ch2=position_peak(peak,[round(x),round(y)],det_im.emccd_im_ch2);
end

%Updated total number of detected photons
sm{tot_phot_det_ch1_idx}=sm{tot_phot_det_ch1_idx}+n_phot_det_ch1;
sm{tot_phot_det_ch2_idx}=sm{tot_phot_det_ch2_idx}+n_phot_det_ch2;







