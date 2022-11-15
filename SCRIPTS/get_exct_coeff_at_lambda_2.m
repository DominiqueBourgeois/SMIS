function updated_spectral_data=get_exct_coeff_at_lambda_2(fluorophore,n_lasers, lasers)

% NAME:
%   get_exct_coeff_at_lambda_2
%
% PURPOSE:
%   This program calculates the exctinction coefficient at particular
%   wavelength from the knowledge of a full spectrum and a reference
%   epsilon at which the spectrum is normalized to 1
%
% INPUTS:
%   fluorophore: the fluorophore
%   lasers: the defined lasers
%
% OUTPUTS:
%	updated_fluorophore: the updated fluorophore
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2012.
%	D.Bourgeois, March 2015.
%	D.Bourgeois, May 2019.
%	D.Bourgeois, March 2021.
%	D.Bourgeois, March 2022. Remove multiple wavelengths possible for a single laser
%

spectral_data=fluorophore.spectral_data; 
updated_spectral_data=spectral_data; % Make a copy of fluorophore
% Define the eps arrays for all lasers during the whole data set
for k=1:spectral_data.n_fluorescent_states
    updated_spectral_data.exc_spectra(k).eps_lasers=zeros(1,n_lasers); 
end
for k=1:spectral_data.n_photoactive_dark_states
    updated_spectral_data.dark_spectra(k).eps_lasers=zeros(1,n_lasers);
end

% Run through all lasers
for i=1:n_lasers
    if lasers(i).power > 0
        for k=1:spectral_data.n_fluorescent_states
            updated_spectral_data.exc_spectra(k).eps_lasers(i)=get_exct_coeff_at_lambda(spectral_data.exc_spectra(k).s,lasers(i).wavelength, updated_spectral_data.exc_spectra(k).eps);
            disp(['Laser: ',num2str(i),' Fluorophore: ', fluorophore.fluorophore_name, ': Extinction coeff of fluorescent state: ', spectral_data.exc_spectra(k).name,' at laser wavelength ',num2str(lasers(i).wavelength),' nm: ', num2str(updated_spectral_data.exc_spectra(k).eps_lasers(i))]);
        end
        for k=1:fluorophore.n_photoactive_dark_states
            updated_spectral_data.dark_spectra(k).eps_lasers(i)=get_exct_coeff_at_lambda(spectral_data.dark_spectra(k).s,lasers(i).wavelength, updated_spectral_data.dark_spectra(k).eps);
            disp(['Laser: ',num2str(i),' Fluorophore: ', fluorophore.fluorophore_name, ': Extinction coeff of dark state: ', spectral_data.dark_spectra(k).name,' at laser wavelength ',num2str(lasers(i).wavelength),' nm: ', num2str(updated_spectral_data.dark_spectra(k).eps_lasers(i))]);
        end
    end
end
