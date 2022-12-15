function c=get_tirf_hilo_correction_factor(z, laser, im_par)

%
% PURPOSE:
%   Get TIRF or HILO correction factor depending on molecule position, laser tirf
%   angle, laser wavelength
%
% INPUTS:
%   z: z position of molecule in raster units
%   laser: laser to analyse
%	im_par: the imaging parameters
%
% OUTPUTS:
%   c: the TIRF amplification factor (c=1 in HILO mode)
%
% MODIFICATION HISTORY:
%	D.Bourgeois, November 2020.

theta=laser.tirf_angle*pi/180; % theta angle in rad

if theta>=im_par.obj.critical_angle   % TIRF mode
    %Formula from Axelrod Meth Enz 2003 eq 4
    c=laser.tirf_amplification*exp(-z*im_par.raster/laser.d);    
else % HILO mode
    theta_sample=asin(im_par.obj.immersion_indice/im_par.obj.sample_indice*sin(theta));
    d=1000*laser.fwhm/tan(theta_sample); % in [nm]
    if z*im_par.raster < d % sm is in the beam
        c=1;
    else % sm is out of the beam, we define as diffraction limited gaussian decay 
        dl=0.61*laser.wavelength/im_par.obj.na; % Size of diffraction limit [nm]
        c=exp(-(z*im_par.raster-d)^2/(2*(dl/2.35)^2));
    end
    
end

