function     im_par=get_camera_signal(im_par)

%
% PURPOSE:
%   Transform number of absorbed photons on detector into an ADU image 
%
% INPUTS:
%	im_par: the imaging parameters
%
% OUTPUTS:
%	im_par: the imaging parameters updated for detector images
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.

% the applied formula is the following
% if N is the number of photons absorbed by the detector, then the measured
% signal in ADU is obtained through:
% Nel_in=N*im_par.det.e2p + *im_par.det.c; % conversion from photons to electrons + adding of spurious charge/clock induced charge
% Nel_out=gamrnd(Nel_in, im_par.det.EMCCD_gain)+ normrnd(0,im_par.det.readout_noise): this is the number of electrons read during the readout process
% ADU=min(floor(Nel_out/im_par.det.e2ADU)+im_par.det.offset, im_par.det.dynamic_range) % # of EMCCD counts
im_par.emccd_im_ch1=min(floor((gamrnd(im_par.emccd_im_ch1*im_par.det.e2p+im_par.det.c, im_par.det.EMCCD_gain)+ normrnd(zeros(size(im_par.emccd_im_ch1)),im_par.det.readout_noise))/im_par.det.e2ADU)+im_par.det.offset, im_par.det.dynamic_range); % # of EMCCD counts
if im_par.two_channel==1
    im_par.emccd_im_ch2=min(floor((gamrnd(im_par.emccd_im_ch2*im_par.det.e2p+im_par.det.c, im_par.det.EMCCD_gain)+ normrnd(zeros(size(im_par.emccd_im_ch1)),im_par.det.readout_noise))/im_par.det.e2ADU)+im_par.det.offset, im_par.det.dynamic_range); % # of EMCCD counts
end



