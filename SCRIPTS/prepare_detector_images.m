function im_par=prepare_detector_images(im_par)
%
% PURPOSE:
%   Prepare empty detector images of the right size
%
% INPUTS:
%	im_par: the imaging parameters
%
% OUTPUTS:
%   im_par: updated im_par for [emccd_im_ch1, emccd_im_ch2, emccd_im_ch1_dl, emccd_im_ch2_dl]
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.

im_par.emccd_im_ch1=double(zeros(im_par.n,im_par.m)); % camera image in channel 1
im_par.emccd_im_ch1_dl=double(zeros(im_par.n,im_par.m)); % diffraction limited image in channel 1 with all SMs

if im_par.two_channel==1
        im_par.emccd_im_ch2=double(zeros(im_par.n,im_par.m)); % camera image in channel 2
        im_par.emccd_im_ch2_dl=double(zeros(im_par.n,im_par.m)); % diffraction limited image in channel 2 with all SMs
    if im_par.single_CCD==1
        im_par.emccd_im_ch12=zeros(im_par.n,2*im_par.m); % Define a dual image for both channels
        im_par.emccd_im_ch12_dl=zeros(im_par.n,2*im_par.m); % Define a dual diffraction limited image for both channels
    end
else
    im_par.emccd_im_ch2=[];
    im_par.emccd_im_ch12_dl=[];
end
