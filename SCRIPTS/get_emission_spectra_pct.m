function sm=get_emission_spectra_pct(sm, sm_par, im_par)

%
% PURPOSE:
%   Get emission spectra and detected number of photons per single
%   molecule
%
% INPUTS:
%   sm: the single molecule
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   sms: the single molecules updated for emission spectra and detected
%   photons
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2022, optimization of parallel computing

%generate the emission spectrum for the SM

%Here are the indices in sm cell array
n_fields=7;
field_idx=num2cell(1:n_fields);
[n_em_idx, n_phot_ch1_idx, n_phot_ch2_idx, tot_n_phot_ch1_idx, ...
    tot_n_phot_ch2_idx, bleached_idx, em_spectrum_idx]=field_idx{:};


%if the molecule was bleached before: nothing to do
if sm{bleached_idx}~=0 && sm{bleached_idx}<im_par.current_frame
    % reset emission spectrum and emitted number of photons if
    % molecules just bleached in preceding frame
    if sm{bleached_idx}==im_par.current_frame-1
        sm{n_phot_ch1_idx}(:)=0;
        if im_par.two_channel==1 ; sm{n_phot_ch2_idx}(:)=0; end
    end
else
    % reset emission spectrum and emitted number of photons
    merged_sm_em_spectrum=[]; % Used if stochastic_spectrum_off=0
    sm{n_phot_ch1_idx}(:)=0;
    if im_par.two_channel==1 ; sm{n_phot_ch2_idx}(:)=0; end

    for i=1:sm_par.n_fluorescent_states
        % only proceed if some photons were emitted
        if sm{n_em_idx}(i)>0
            em_spectrum=sm_par.spectral_data.em_spectra(i).s; % Reference emission spectrum

            if ~im_par.stochastic_spectrum_off
                %Initialize .em_spectrum field if necessary
                if isempty(sm{em_spectrum_idx})
                    sm{em_spectrum_idx}.frames=im_par.current_frame;
                    sm{em_spectrum_idx}.spectra=[];
                end

                %Get the spectrum
                sm_em_spectrum = get_sm_emitted_spectrum(sm{n_em_idx}(i), em_spectrum, im_par);

                %Add the frame tag
                if sm{em_spectrum_idx}.frames(end)<im_par.current_frame
                    sm{em_spectrum_idx}.frames=[sm{em_spectrum_idx}.frames, im_par.current_frame];
                end

                %Combine the spectra from all fluorescent states
                merged_sm_em_spectrum=merge_em_spectrum(merged_sm_em_spectrum,sm_em_spectrum);


            else
                sm_em_spectrum = [em_spectrum(:,1), sm{n_em_idx}(i)*em_spectrum(:,2)/sum(em_spectrum(:,2))];
            end

            %filter that spectrum towards channels 1 & 2
            [sm{n_phot_ch1_idx}(i), sm{n_phot_ch2_idx}(i)] = filter_em_spectrum(sm_em_spectrum, sm_par.filter_profiles(i), im_par);

            %Update total # of photons
            sm{tot_n_phot_ch1_idx}(i)=sm{tot_n_phot_ch1_idx}(i)+sm{n_phot_ch1_idx}(i);
            sm{tot_n_phot_ch2_idx}(i)=sm{tot_n_phot_ch2_idx}(i)+sm{n_phot_ch2_idx}(i);
        end
    end

    if ~im_par.stochastic_spectrum_off && ~isempty(merged_sm_em_spectrum)
        sm{em_spectrum_idx}.spectra=[sm{em_spectrum_idx}.spectra,merged_sm_em_spectrum];
    end
end
