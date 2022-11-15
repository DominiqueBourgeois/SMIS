function updated_spectral_data=get_exct_coeff_at_lambda_3(fluorophore,n_lasers, lasers)

% NAME:
%   get_exct_coeff_at_lambda_3
%
% PURPOSE:
%   This program calculates the exctinction coefficient at particular
%   wavelength from the knowledge of a full spectrum and a reference
%   epsilon at which the spectrum is normalized to 1 in a pH dependant
%   manner
%
% INPUTS:
%   dye: the fluorophore
%   lasers: the defined lasers
%
% OUTPUTS:
%	updated_spectral_data: the updated fluorophore
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2020.
%	D.Bourgeois, March 2022. Remove multiple wavelengths possible for a single laser
%
spectral_data=fluorophore.spectral_data; 
photoactive_dark_states=fluorophore.photoactive_dark_states;

updated_spectral_data=spectral_data; % Make a copy of dye spectral data
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
            % Correct for pH dependant reference eps
            updated_spectral_data.exc_spectra(k).eps=fluorophore.fluorescent_fraction(k)*spectral_data.exc_spectra(k).eps;
            updated_spectral_data.exc_spectra(k).eps_lasers(i)=get_exct_coeff_at_lambda(spectral_data.exc_spectra(k).s,lasers(i).wavelength, updated_spectral_data.exc_spectra(k).eps);
            disp(['Laser: ',num2str(i),' Dye: ', fluorophore.fluorophore_name, ': Extinction coeff of fluorescent state: ', spectral_data.exc_spectra(k).name,' at laser wavelength ',num2str(lasers(i).wavelength),' nm: ', num2str(updated_spectral_data.exc_spectra(k).eps_lasers(i))]);
        end
        for k=1:spectral_data.n_photoactive_dark_states
            %check if the dark state is associated to a flurescent state
            %with pH equilibrium           
            w_assoc_to_fluo=find(photoactive_dark_states(k)==fluorophore.associated_dark_states,1);
            if ~isempty(w_assoc_to_fluo)
                fr_dark=1-fluorophore.fluorescent_fraction(w_assoc_to_fluo);                
                updated_spectral_data.dark_spectra(k).eps=fr_dark*spectral_data.dark_spectra(k).eps;
                updated_spectral_data.dark_spectra(k).eps_lasers(i)=get_exct_coeff_at_lambda(spectral_data.dark_spectra(k).s,lasers(i).wavelength, updated_spectral_data.dark_spectra(k).eps);
                disp(['Laser: ',num2str(i),' Dye: ', fluorophore.fluorophore_name, ': Extinction coeff of neutral state: ', spectral_data.dark_spectra(k).name,' at laser wavelength ',num2str(lasers(i).wavelength),' nm: ', num2str(updated_spectral_data.dark_spectra(k).eps_lasers(i))]);
            else
                updated_spectral_data.dark_spectra(k).eps_lasers(i)=get_exct_coeff_at_lambda(spectral_data.dark_spectra(k).s,lasers(i).wavelength, spectral_data.dark_spectra(k).eps);
                disp(['Laser: ',num2str(i),' Dye: ', fluorophore.fluorophore_name, ': Extinction coeff of dark state: ', spectral_data.dark_spectra(k).name,' at laser wavelength ',num2str(lasers(i).wavelength),' nm: ', num2str(updated_spectral_data.dark_spectra(k).eps_lasers(i))]);
            end
        end
    end
end