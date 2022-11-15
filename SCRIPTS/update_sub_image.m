function [sm, im_par] = update_sub_image(x,y,z,sm, im_par, sm_par, channel)
%
% NAME:
%	update_sub_image
% PURPOSE:
%   Case of a diffuse PSF: the PSF of a single molecule is approximated by
%   the superposition of several PSF's from the individual sub_xyz
%   coordinates
% INPUTS:
%   x, y, z: Multiple coordinates of the moving single molecule on detector image [raster
%   units]
%   sm: the single molecule (with coordinates on high-resolution image) in raster units
%	im_par: the imaging parameters
%	sm_par: the sm parameters
%   channel: set to channel number (1 for channel 1, etc)
%
% OUTPUTS:
%   im_par: the imaging parameters updated for detector image with the contribution from the new
%   Moving mocalization
%   sm: The updated single molecule
%
% MODIFICATION HISTORY:
%	D.Bourgeois, March 2012. Modified September 2019
%	D.Bourgeois, May 2021. Introduce sm.n_phot_det_ch1 and sm.n_phot_det_ch2;
%	D.Bourgeois, July 2021. Keep track of sm.n_phot_det_ch1 and sm.n_phot_det_ch2 over whole data set;

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
    if channel==1 && sm.n_phot_ch1(i)>0
        %get the number of photons actually reaching the detector
        n_phot_det_ch1=n_phot_det_ch1+get_emitted_photons(sm.n_phot_ch1(i), sm.theta, im_par, sm_par);
    elseif channel==2 && sm.n_phot_ch2(i)>0 %same for channel 2
        %get the number of photons actually reaching the detector
        n_phot_det_ch2=n_phot_det_ch2+get_emitted_photons(sm.n_phot_ch2(i), sm.theta, im_par, sm_par);
    end
end

% treat channel 1 first
if channel==1 && n_phot_det_ch1>0
    %update sm.n_phot_det_ch1
    sm.n_phot_det_ch1=[sm.n_phot_det_ch1,n_phot_det_ch1];
    sm.frames_on_ch1=[sm.frames_on_ch1,im_par.current_frame];

    for k=1:numel(sm.sub_x) % Go for all sub-positions
        if sm.fr_t_on(i,k)>0 % Only calculate the PSF if there are some photons
            if im_par.simul_3D==0
                peak=get_2Dpeak(x(k),y(k),round(n_phot_det_ch1*sm.fr_t_on(i,k)) ,sm_par,im_par,i,channel);
            else
                peak=get_3Dpeak(x(k),y(k),z(k), round(n_phot_det_ch1*sm.fr_t_on(i,k)), sm_par,im_par,i,channel);
            end
            %Position this spot into the whole image
            im_par.emccd_im_ch1=position_peak(peak,[round(x(k)),round(y(k))],im_par.emccd_im_ch1);
        end
    end
elseif channel==2 && n_phot_det_ch2>0 %same for channel 2
    %update sm.n_phot_det_ch2
    sm.n_phot_det_ch2=[sm.n_phot_det_ch2,n_phot_det_ch2];
    sm.frames_on_ch2=[sm.frames_on_ch2,im_par.current_frame];

    for k=1:numel(sm.sub_x) % Go for all sub-positions
        if sm.fr_t_on(i,k)>0 % Only calculate the PSF if there are some photons
            if im_par.simul_3D==0
                peak=get_2Dpeak(x(k),y(k),round(n_phot_det_ch2*sm.fr_t_on(i,k)),sm_par,im_par,i,channel);
            else
                peak=get_3Dpeak(x(k),y(k),z(k),round(n_phot_det_ch2*sm.fr_t_on(i,k)),sm_par,im_par,i,channel);
            end
            %Position this spot into the whole image
            im_par.emccd_im_ch2=position_peak(peak,[round(x(k)),round(y(k))],im_par.emccd_im_ch2);
        end
    end
end

%Updated total number of detected photons
sm.tot_phot_det_ch1=sm.tot_phot_det_ch1+n_phot_det_ch1;
sm.tot_phot_det_ch2=sm.tot_phot_det_ch2+n_phot_det_ch2;





