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
%   par: a structure containing some imaging parameters (n, m, raster).
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


c=[par.m/2+0.5 par.n/2+0.5];
[xpix,ypix] = meshgrid(1:par.m,1:par.n);
xpix=xpix-c(1);
ypix=ypix-c(2);

% run over all lasers
for i=1:n_lasers
    if lasers(i).power>0
        switch lasers(i).mode
            case char('Gaussian')
                sigmax = 1000*lasers(i).fwhm/par.raster/2.35; % in pixels
                sigmay = 1000*lasers(i).fwhm/par.raster/2.35; % in pixels
                u=((xpix).*(xpix))/(2*sigmax^2)+((ypix).*(ypix))/(2*sigmay^2);
                beam_profile=exp(-u);
                %normalize to the actual total energy [J]
                %maybe the laser spot extends far out of the image, so we need to normalize
                %relative to the full laser spot, which we take at 10 sigmas
                c_ref=[round(10*sigmax)/2+0.5 round(10*sigmay)/2+0.5];
                [xpix_ref,ypix_ref] = meshgrid(1:round(10*sigmax),1:round(10*sigmay));
                xpix_ref=xpix_ref-c_ref(1);
                ypix_ref=ypix_ref-c_ref(2);
                u_ref=((xpix_ref).*(xpix_ref))/(2*sigmax^2)+((ypix_ref).*(ypix_ref))/(2*sigmay^2);
                beam_profile_ref=exp(-u_ref);
                beam_profile=beam_profile/sum(sum(beam_profile_ref))*lasers(i).power*1e-03; % Available power [W] in camera FOV
                %Calculate power density [W/cm²] at center
                %raster in nm => 1e+7 to go to cm ; duration in ms => 1e+3 to go to seconds
                lasers(i).power_density=max(max(beam_profile))/((par.raster*1e-7)^2);
                disp(['Power density of laser:', num2str(i),' at center [W/cm^2] (100%): ', num2str(lasers(i).power_density)]);               
            case char('Flat')
                beam_profile=ones(par.n,par.m);
                sxy = 1000*lasers(i).fwhm/par.raster; % laser size in pixels
                if sxy<par.n % If laser size smaller than FOV, set to 0 the missing part
                    beam_profile(1:round((par.n-sxy)/2),:)=0;
                    beam_profile(par.n-round((par.n-sxy)/2):par.n,:)=0;
                end
                if sxy<par.m % If laser size smaller than FOV, set to 0 the missing part
                    beam_profile(:,1:round((par.m-sxy)/2))=0;
                    beam_profile(:,par.m-round((par.m-sxy)/2):par.m)=0;
                end
                
                %normalize to the actual total energy [J]
                beam_profile=beam_profile/(sxy^2)*lasers(i).power*1e-03; % Available power [W] in camera FOV
                %Calculate power density [W/cm²] at center
                %raster in nm => 1e+7 to go to cm ; duration in ms => 1e+3 to go to seconds
                lasers(i).power_density=max(max(beam_profile))/((par.raster*1e-7)^2);
                disp(['Power density of laser:', num2str(i),' [W/cm^2] (100%): ', num2str(lasers(i).power_density)]);
            otherwise
                error('Error, laser shape not recognized ! Can be ''Gaussian'' or ''Flat''!')
        end
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

