function ens=get_ensemble_fluorescence(ens, sm_par, im_par)

%
% PURPOSE:
%   Get the ensemble fluorescence for all fluorophores during current frame
%
% INPUTS:
%   ens: the ensemble data
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   ens: The updated ensemble data
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2021.

n_fluorophores=numel(sm_par);

for i=1:n_fluorophores
    %First we have to filter the signal through the emission filters in the channels
    fr_detected=zeros(2,sm_par(i).n_fluorescent_states);
    for k=1:sm_par(i).n_fluorescent_states
        %Normalize the emission spectrum
        em_spectrum=sm_par(i).spectral_data.em_spectra(k).s;
        em_spectrum = [em_spectrum(:,1), em_spectrum(:,2)/sum(em_spectrum(:,2))];
        
        %Get the fraction of the spectrum that goes through channels 1 & 2

        filters=sm_par(i).filter_profiles(k);
        [fr_detected(1,k), fr_detected(2,k)] = filter_em_spectrum_ensemble(em_spectrum, filters, im_par);
    end
    
    %Second get the signal from the fluorescent states in each
    %subpopulation in each channel
    for k=1:sm_par(i).n_subpop
        F=ens(i).sp(k).det_p(sm_par(i).fluorescent_states,im_par.current_frame); % The state population for the fluorescent states
        % Now we multiply this by number of photons per frame times quantum yield times fluorescent fraction
        S=F'.*(ens(i).sp(k).Ns(end,sm_par(i).fluorescent_states)*im_par.frametime*1e-3).*sm_par(i).quantum_yield.*sm_par(i).fluorescent_fraction;
        ens(i).sp(k).S_ch1(:,im_par.current_frame)=S'.*fr_detected(1,:)';
        if im_par.two_channel==1
            ens(i).sp(k).S_ch2(:,im_par.current_frame)=S'.*fr_detected(2,:)';
        end
    end
    
    %Get the total signal over the different subpopulations
    for k=1:sm_par(i).n_subpop
        ens(i).det_p(:,im_par.current_frame)=ens(i).det_p(:,im_par.current_frame)+ens(i).sp(k).det_p(:,im_par.current_frame);
        ens(i).S_ch1(:,im_par.current_frame)=ens(i).S_ch1(:,im_par.current_frame)+ens(i).sp(k).S_ch1(:,im_par.current_frame);
        if im_par.two_channel==1
            ens(i).S_ch2(:,im_par.current_frame)=ens(i).S_ch2(:,im_par.current_frame)+ens(i).sp(k).S_ch2(:,im_par.current_frame);
        end
    end
    %Get the average
    ens(i).det_p(:,im_par.current_frame)=ens(i).det_p(:,im_par.current_frame)/sm_par(i).n_subpop;
    ens(i).S_ch1(:,im_par.current_frame)=ens(i).S_ch1(:,im_par.current_frame)/sm_par(i).n_subpop;
    if im_par.two_channel==1
        ens(i).S_ch2(:,im_par.current_frame)=ens(i).S_ch2(:,im_par.current_frame)/sm_par(i).n_subpop;
    end
end
