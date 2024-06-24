function epsilon = get_exct_coeff_at_lambda(spectrum, lambda, epsilon_ref)
% NAME:
%   get_exct_coeff_at_lambda
%
% PURPOSE:
%   This program calculates the exctinction coefficient at a particular
%   wavelength from the knowledge of a full spectrum and a reference
%   epsilon at which the spectrum is normalized to 1
%
% INPUTS:
%   spectrum: the spectrum
%   lambda: the wavelength [nm]
%   epsilon_ref: the reference epsilon at which the spectrum is normalized to 1
%
% OUTPUTS:
%	epsilon: the calculated epsilon
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2012.
%	D.Bourgeois, March 2015.
%	D.Bourgeois, October 2020. Handle wavelengths very remote from
%	available spectra by linear interpolation
%
n=length(lambda);
epsilon=zeros(n,1);
min_lambda_discrepancy=2; % [nm] Tolerable distance to lambda for spectrum data 
for k=1:n
    tmp_var=abs(spectrum(:,1)-lambda(k));
    w_lambda=find(tmp_var == min(tmp_var),1); %index of the position where to normalize
    if abs(spectrum(w_lambda,1)-lambda(k))>min_lambda_discrepancy
        %Try to interpolate linearly
        tmp_var=spectrum(:,1)-lambda(k);
        if ~isempty(tmp_var(tmp_var>0))
            w1_lambda=find(tmp_var==min(tmp_var(tmp_var>0)),1); %index of the position for closest wavelength > lambda(k)
        else
            w1_lambda=[];
        end
        if ~isempty(tmp_var(tmp_var<0))
            w2_lambda=find(tmp_var==max(tmp_var(tmp_var<0)),1); %index of the position for closest wavelength < lambda(k)
        else
            w2_lambda=[];
        end

        if ~isempty(w1_lambda) && ~isempty(w2_lambda) % estimate in between        
            sl=(spectrum(w1_lambda,2)-spectrum(w2_lambda,2))/(spectrum(w1_lambda,1)-spectrum(w2_lambda,1)); %slope
            of=(spectrum(w2_lambda,2)*spectrum(w1_lambda,1)-spectrum(w1_lambda,2)*spectrum(w2_lambda,1))/(spectrum(w1_lambda,1)-spectrum(w2_lambda,1)); %offset
            epsilon(k)=epsilon_ref*(sl*lambda(k)+of);
        else % Then impossible to estimate excinction coeff
            SMISMessage1=['Spectrum does not allow estimating exctinction coeff at: ',num2str(lambda(k)),' nm. Check fluorophore definition !'];
            SMISMessage2=['Will use epsilon = 0 at wavelength: ',num2str(lambda(k)), ' nm'];
            SMISMessage=[SMISMessage1 newline SMISMessage2];
            disp(SMISMessage);
            %warndlg(SMISMessage,'Warning')
            %uiwait
            return
        end
    else % Normal case
        epsilon(k)=epsilon_ref*spectrum(w_lambda,2);
    end
end
end