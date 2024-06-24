function im_par=get_tx_background_ch(im_par, tx_pattern, attenuation_factor, lasers, channel)

% PURPOSE:Les
%   Get background on Channel 1 or 2
%
% INPUTS:
%	im_par: the imaging parameters
%   tx_image: the texture binary image
%	attenuation_factor: the attenuation factor due to accumulated dose and
%   lasers: the lasers
%   channel: the channel number
%
% OUTPUTS:
%	im_par: the imaging parameters updated for detector images in the
%	specified channel
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.
%	D.Bourgeois, June 2021: Make background calculation dependent on laser
%	power density and checked that laser is on during frame time
%	D.Bourgeois, September 2021: Make background calculation per s instead
%	of per frame

if channel>2; message('Channel does not exist !'); return; end
if channel==1
    tx_bg_p = im_par.bg.textured_bg_ch1*(im_par.raster*0.01)^2; %background per pixel
    tx_bg_laser_sensitivity=im_par.bg.textured_bg_laser_sensitivity_ch1;
else
    tx_bg_p = im_par.bg.textured_bg_ch2*(im_par.raster*0.01)^2; %background per pixel
    tx_bg_laser_sensitivity=im_par.bg.textured_bg_laser_sensitivity_ch2;
end

n_lasers=size(lasers,2); % # of lasers

if tx_bg_p>0
    if im_par.bg.adjust_textured_bg==1 
        % calculate the background fluorescence as:
        % bg_fluo=attenuation_factor*intensity
        % the intensity is calculated in such a way that the total number
        % of detected photons in the first frame for all lasers is bg_p
        % photons per pixel (in the central pixel)in channel 1.
        tx_bg_fluo=zeros(im_par.n,im_par.m);
        
        for k=1:n_lasers % add contributions from all lasers (in # of photons/pixel)
            % Evaluate if we are in TIRF or HILO mode and reduce background
            % accordingly
            if lasers(k).tirf==1
                theta=lasers(k).tirf_angle*pi/180; % theta angle in rad
                sample_thickness=im_par.nz*im_par.raster; % [nm]
                if theta>=im_par.obj.critical_angle   % TIRF mode
                    %get integrated laser power density over the TIRF
                    %thickness and normalize to the sample thickness
                    d=lasers(k).d;
                    tirf_eff_sample_thickness=lasers(k).tirf_amplification*d*(1-exp(-sample_thickness/d));
                    tirf_hilo_corr=min([1,tirf_eff_sample_thickness/sample_thickness]); % If TIRF beam is thicker than sample, correction factor should be 1                                  
                else % HILO mode
                    theta_sample=asin(im_par.obj.immersion_indice/im_par.obj.sample_indice*sin(theta));
                    hilo_eff_sample_thickness=1000*lasers(k).fwhm/tan(theta_sample); % in [nm]
                    tirf_hilo_corr=min([1,hilo_eff_sample_thickness/sample_thickness]); % If HILO beam is thicker than sample, correction factor should be 1                 
                end
            else
                tirf_hilo_corr=1; % No correction to apply
            end
                        
            tx_bg_fluo=tx_bg_fluo+tirf_hilo_corr*tx_bg_laser_sensitivity(k)*...
                lasers(k).on_during_frametime*...
                lasers(k).power_density/100*...
                lasers(k).sequence(im_par.current_frame)/100*...
                lasers(k).beam_profile/lasers(k).max_beam_profile*...
                im_par.frametime*1e-3;
            
        end
        tx_bg_fluo=tx_bg_fluo.*tx_pattern; % here multiply by the pattern
        % the attenuated signal is then
        tx_bg_fluo=imnoise(uint16(tx_bg_p*attenuation_factor.*tx_bg_fluo),'poisson');
        %Note that there are rounding errors apparently that makes the
        %noise pattern not very smooth along the beam profile
    else % Simple estimation of background
        tx_bg_fluo=zeros(im_par.n,im_par.m);
        
        for k=1:n_lasers % add contributions from all lasers (in # of photons/pixel)
            
            tx_bg_fluo=tx_bg_fluo+1/n_lasers*...
                lasers(k).on_during_frametime*... % Only if laser present during frametime
                lasers(k).power_density/100*... % [W/cm²]Power
                im_par.frametime*1e-3; % duration [s]
        end
        % the attenuated signal is then
        tx_bg_fluo=uint16(tx_bg_p.*tx_bg_fluo.*tx_pattern);
        tx_bg_fluo=imnoise(tx_bg_fluo,'poisson');
        
        %tx_bg_fluo=imnoise(uint16(tx_bg_p*tx_pattern),'poisson');
    end      
end    
    %%
    
    
    if channel==1
        im_par.emccd_im_ch1=im_par.emccd_im_ch1+double(tx_bg_fluo);
    else
        im_par.emccd_im_ch2=im_par.emccd_im_ch2+double(tx_bg_fluo);
    end
end

