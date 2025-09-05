function lasers = get_beam_profile(n_lasers, lasers, par)
% NAME:
%	GET_BEAM_PROFILE
%
% PURPOSE:
%	Derive a Gaussian beam profile
%
% CATEGORY:
%	Signal, image processing.
%
% INPUTS:
%   n_lasers: number of lasers
%	laser: a structure containing laser parameters (power, duration, fwhm, wavelength).
%   par: a structure containing some imaging parameters (raster).
%
% OUTPUTS:
%	laser: updated laser structure with beam_profile, power density, max power density
%   laser.beam_profile in units of [# photons/second]
%
% COMMON BLOCKS:
%	None.
%
% SIDE EFFECTS:
%	None.
%
% RESTRICTIONS:
%	None.
%
% MODIFICATION HISTORY:
%	D.Bourgeois, April 2011.
%	D.Bourgeois, April 2013 added return of max_beam_profile
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15
%	D.Bourgeois, June 2020: beam_profile was calculated at first frame. Was
%	inconsistent with get_number_of_absorbed_photons
%	D.Bourgeois, November 2020: option for flat beam_profile 
%	D.Bourgeois, June 2022: remove global variables

%Some variables
planck = 6.62e-34; % [J.s]
speed_of_light = 299790000; %[m/s]

% run over all lasers
for i=1:n_lasers
    if lasers(i).power>0
        beam_profile=lasers(i).beam_profile;
        %Calculate power density [W/cm²] at max of profile
        %raster in nm => 1e+7 to go to cm ; duration in ms => 1e+3 to go to seconds
        lasers(i).power_density=max(max(beam_profile))/((par.raster*1e-7)^2);
        disp(['Power density of laser at peak location:', num2str(i),' [W/cm^2] (100%): ', num2str(lasers(i).power_density)]);

        %Convert from W to number of photons/s.
        %Energy of one photon [J/ph]
        e_phot=planck*speed_of_light/(lasers(i).wavelength*1e-09);
        lasers(i).beam_profile=beam_profile/e_phot;
        
        %Get also the maximum value # of photons/s in central pixel
        lasers(i).max_beam_profile=max(max(lasers(i).beam_profile));
        
        %Get also the maximum energy J/s in central pixel
        lasers(i).max_beam_profile_watt=lasers(i).max_beam_profile*e_phot;
    end
end

end

