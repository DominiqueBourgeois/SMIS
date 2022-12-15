function     det_im=get_camera_signal_pct(det_im,im_par)

%
% PURPOSE:
%   Transform number of absorbed photons on detector into an ADU image 
%
% INPUTS:
%	det_im: the detector images 
%	im_par: the imaging parameters
%
% OUTPUTS:
%	det_im: the updated detector images 
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.
%	D.Bourgeois, September 2022. Decouple det_im from im_par to reduce size
%	of im_par as a broadcast variable in pct


% the applied formula is the following
% if N is the number of photons absorbed by the detector, then the measured
% signal in ADU is obtained through:
% Nel_in=N*im_par.det.e2p + *im_par.det.c; % conversion from photons to electrons + adding of spurious charge/clock induced charge
% Nel_out=gamrnd(Nel_in, im_par.det.EMCCD_gain)+ normrnd(0,im_par.det.readout_noise): this is the number of electrons read during the readout process
% ADU=min(floor(Nel_out/im_par.det.e2ADU)+im_par.det.offset, im_par.det.dynamic_range) % # of EMCCD counts
det_im.emccd_im_ch1=min(floor((gamrnd(det_im.emccd_im_ch1*im_par.det.e2p+im_par.det.c, im_par.det.EMCCD_gain)+ normrnd(zeros(size(det_im.emccd_im_ch1)),im_par.det.readout_noise))/im_par.det.e2ADU)+im_par.det.offset, im_par.det.dynamic_range); % # of EMCCD counts
if im_par.two_channel==1
    det_im.emccd_im_ch2=min(floor((gamrnd(det_im.emccd_im_ch2*im_par.det.e2p+im_par.det.c, im_par.det.EMCCD_gain)+ normrnd(zeros(size(det_im.emccd_im_ch1)),im_par.det.readout_noise))/im_par.det.e2ADU)+im_par.det.offset, im_par.det.dynamic_range); % # of EMCCD counts
end



