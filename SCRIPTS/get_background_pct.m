function     [det_im, lasers]=get_background_pct(det_im, im_par, lasers, tx_pattern)

%
% PURPOSE:Les
%   Add background on detector images
%
% INPUTS:
%	det_im: the detector images 
%	im_par: the imaging parameters
%	lasers: the lasers
%   tx_pattern: 2D image of textured background. [] if no textured background used
%
% OUTPUTS:
%	det_im: the updated detector images 
%	lasers: the lasers Updated for Accumulated doses
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.
%	D.Bourgeois, September 2022, introduce det_im
%	D.Bourgeois, February 2023, use tx_pattern as new input variable to avoid rereading at every frame

% if there is some attenuation, calculate the accumulated dose to the sample
if (im_par.bg.adjust_bg==1 && im_par.bg.bg_decay_rate>0) || ...
        (im_par.bg.adjust_textured_bg==1 && im_par.bg.textured_bg_decay_rate>0)
    [dose, lasers]=update_lasers_accumulated_dose(lasers, im_par);
end

if im_par.bg.adjust_bg==1 && im_par.bg.bg_decay_rate>0
    % The decay rate per pixel is
    rate_p=im_par.bg.bg_decay_rate*(im_par.raster*0.01)^2;
    attenuation_factor=exp(-rate_p*dose);
else
     attenuation_factor=1; % here a single constant is enough
end

if im_par.bg.adjust_textured_bg==1 && im_par.bg.textured_bg_decay_rate>0
    % The decay rate per pixel is
    tx_rate_p=im_par.bg.textured_bg_decay_rate*(im_par.raster*0.01)^2;
    tx_attenuation_factor=exp(-tx_rate_p*dose);
else
    tx_attenuation_factor=1; % here a single constant is enough
end

%Channel 1
if im_par.bg.bg_ch1>0
    det_im=get_background_ch_pct(det_im, im_par, attenuation_factor, lasers, 1);
end
if im_par.bg.add_textured_bg==1 && im_par.bg.textured_bg_ch1>0
    det_im=get_tx_background_ch_pct(det_im, im_par, tx_pattern, tx_attenuation_factor, lasers, 1);
end

%Channel 2
if im_par.two_channel==1 && im_par.bg.bg_ch2>0
    det_im=get_background_ch_pct(det_im, im_par, attenuation_factor, lasers, 2);
end
if im_par.bg.add_textured_bg==1 && im_par.two_channel==1 && im_par.bg.textured_bg_ch2>0
    det_im=get_tx_background_ch_pct(det_im, im_par, tx_pattern, tx_attenuation_factor, lasers, 2);
end

