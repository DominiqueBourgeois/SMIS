function sms=get_emission_spectra(sms, sm_par, im_par)

%
% PURPOSE:
%   Get emission spectra and detected number of photons per single
%   molecule
%
% INPUTS:
%   sms: the single molecules (with coordinates on high-resolution image) in raster units
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   sms: the single molecules updated for emission spectra and detected
%   photons
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.

%generate the emission spectrum for the SM

for k=1:sm_par.n_mol_eff
    %if the molecule was bleached before: nothing to do
    if sms.sm(k).bleached~=0 && sms.sm(k).bleached<im_par.current_frame
        % reset emission spectrum and emitted number of photons if
        % molecules just bleached in preceding frame
        if sms.sm(k).bleached==im_par.current_frame-1
            sms.sm(k).n_phot_ch1(:)=0;
            if im_par.two_channel==1 ; sms.sm(k).n_phot_ch2(:)=0; end
        end
    else
        % reset emission spectrum and emitted number of photons
        merged_sm_em_spectrum=[]; % Used if stochastic_spectrum_off=0
        sms.sm(k).n_phot_ch1(:)=0;
        if im_par.two_channel==1 ; sms.sm(k).n_phot_ch2(:)=0; end
        
        for i=1:sm_par.n_fluorescent_states
            % only proceed if some photons were emitted
            if sms.sm(k).n_em(i)>0
                em_spectrum=sm_par.spectral_data.em_spectra(i).s; % Reference emission spectrum

                if ~im_par.stochastic_spectrum_off
                    %Initialize .em_spectrum field if necessary
                    if isempty(sms.sm(k).em_spectrum)
                        sms.sm(k).em_spectrum.frames=im_par.current_frame;
                        sms.sm(k).em_spectrum.spectra=[];
                    end

                    %Get the spectrum
                    sm_em_spectrum = get_sm_emitted_spectrum(sms.sm(k).n_em(i), em_spectrum, im_par);
                    
                    %Add the frame tag 
                    if sms.sm(k).em_spectrum.frames(end)<im_par.current_frame
                        sms.sm(k).em_spectrum.frames=[sms.sm(k).em_spectrum.frames, im_par.current_frame];
                    end
                    
                    %Combine the spectra from all fluorescent states
                    merged_sm_em_spectrum=merge_em_spectrum(merged_sm_em_spectrum,sm_em_spectrum);


                else
                    sm_em_spectrum = [em_spectrum(:,1), sms.sm(k).n_em(i)*em_spectrum(:,2)/sum(em_spectrum(:,2))];
                end
                
                %filter that spectrum towards channels 1 & 2
                [sms.sm(k).n_phot_ch1(i), sms.sm(k).n_phot_ch2(i)] = filter_em_spectrum(sm_em_spectrum, im_par);
                
                %Update total # of photons
                sms.sm(k).tot_n_phot_ch1(i)=sms.sm(k).tot_n_phot_ch1(i)+sms.sm(k).n_phot_ch1(i);
                sms.sm(k).tot_n_phot_ch2(i)=sms.sm(k).tot_n_phot_ch2(i)+sms.sm(k).n_phot_ch2(i);
            end
        end

        if ~im_par.stochastic_spectrum_off && ~isempty(merged_sm_em_spectrum)
            sms.sm(k).em_spectrum.spectra=[sms.sm(k).em_spectrum.spectra,merged_sm_em_spectrum];
        end
    end
end
