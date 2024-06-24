function det_im=prepare_detector_images_pct(im_par)
%
% PURPOSE:
%   Prepare empty detector images of the right size
%
% INPUTS:
%	im_par: the imaging parameters
%
% OUTPUTS:
%   det_im: updated structure containing [emccd_im_ch1, emccd_im_ch2, emccd_im_ch1_dl, emccd_im_ch2_dl]
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.
%	D.Bourgeois, September 2022. Decouple det_im from im_par to reduce size
%	of im_par as a broadcast variable in pct


det_im.emccd_im_ch1=double(zeros(im_par.n,im_par.m)); % camera image in channel 1
det_im.emccd_im_ch1_dl=double(zeros(im_par.n,im_par.m)); % diffraction limited image in channel 1 with all SMs

if im_par.two_channel==1
        det_im.emccd_im_ch2=double(zeros(im_par.n,im_par.m)); % camera image in channel 2
        det_im.emccd_im_ch2_dl=double(zeros(im_par.n,im_par.m)); % diffraction limited image in channel 2 with all SMs
    if im_par.single_CCD==1
        det_im.emccd_im_ch12=zeros(im_par.n,2*im_par.m); % Define a dual image for both channels
        det_im.emccd_im_ch12_dl=zeros(im_par.n,2*im_par.m); % Define a dual diffraction limited image for both channels
    end
else
    det_im.emccd_im_ch2=[];
    det_im.emccd_im_ch12_dl=[];
end
