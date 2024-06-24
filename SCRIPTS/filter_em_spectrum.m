function [n1, n2] = filter_em_spectrum(spectrum, filter_profiles, im_par)
% NAME:
%	filter_em_spectrum
%
% PURPOSE:
%	Filter spectrum by 2 filters for 2-channel experiments
%
% CATEGORY:
%	Signal, image processing.
%
% CALLING SEQUENCE:
%	[n1, n2] = filter_em_spectrum(spectrum, par)
%
% INPUTS:
%   spectrum: the spectrum to be filtered
%	filter_profiles: the filter_profiles defined in get_sm_filter_profiles
%	im_par: SMIS imaging parameters
%
% OUTPUTS:
%	n1, n2: Number of photons emitted on each channel
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
%	D.Bourgeois, April 2011, September 2017, September 2019.
%	D.Bourgeois, September 2022 Apply Poisonian noise after going through
%	the filters (except if im_par.stochastic_spectrum_off=0). If Poissonian noise would be applied before (at the level
%	of the # of emitted photons), then the ratio between photons reaching
%	channel 1 and 2 would not be affected by Poissonian noise, which would
%	be incorrect for eg spectral demixing approaches.

% l=spectrum(:,1); % The wavelengths
s=spectrum(:,2); % The spectral values

% f1 = im_par.filters.ch1_filter; % Get filter 1
% n_bands1=size(f1,1); % # of bands in the filter
% n1=0;

if im_par.two_channel==0
%         for i=1:n_bands1
%             if f1(i,1)<f1(i,2) % First filter definition mode [lambda_min, lambda_max, transmission]
%                 w_ch1 = find((l >= f1(i,1) & l <= f1(i,2)));
%             else % Second filter definition mode [lambda_central, bandpass, transmission]
%                 w_ch1 = find((l >= (f1(i,1)-f1(i,2)/2)) & (l <= (f1(i,1)+f1(i,2)/2)));
%             end
%             if ~isempty(w_ch1)
%                 n1 = n1+fix(sum(s(w_ch1))*f1(i,3));
%             end
%         end

    n1=sum(s.*filter_profiles.ch1);

    %Apply Poissonian law if needed
    if im_par.stochastic_spectrum_off==1
        n1=poissrnd(n1);
    end

    n2=0;


else % two channel experiment
    %     f2 = im_par.filters.ch2_filter; % Get filter 2
    %     n_bands2=size(f2,1); % # of bands in the filter

    % First consider eventual dichroic
    if im_par.filters.add_dichroic==1
        %         df=im_par.filters.dichroic_filter;
        %         if min(l)<(df(1)-df(2)/2) && max(l)>(df(1)+df(2)/2) % General case
        %             xdf=[min(l),df(1)-df(2)/2,df(1)+df(2)/2, max(l)];
        %             ydf=[df(3), df(3), df(4),df(4)];
        %         elseif min(l)<(df(1)-df(2)/2) && max(l)<=(df(1)+df(2)/2)
        %             xdf=[min(l),df(1)-df(2)/2,df(1)+df(2)/2];
        %             ydf=[df(3), df(3), df(4)];
        %         elseif min(l)>=(df(1)-df(2)/2) && max(l)>(df(1)+df(2)/2)
        %             xdf=[df(1)-df(2)/2,df(1)+df(2)/2,max(l)];
        %             ydf=[df(3), df(4), df(4)];
        %         else
        %             xdf=[df(1)-df(2)/2,df(1)+df(2)/2];
        %             ydf=[df(3), df(4)];
        %         end
        %         f=interp1(xdf,ydf,l, "linear");

        s2=s.*filter_profiles.df; % spectrum reaching channel 2
        s=s.*(1-filter_profiles.df); % spectrum reaching channel 1
        
        %s2=[s(l<=df(1))*df(2); s(l>df(1))*df(3)]; % spectrum reaching channel 2
        %s=[s(l<=df(1))*(1-df(2)); s(l>df(1))*(1-df(3))]; % spectrum reaching channel 1
    else
        s2=s;
    end

    %First channel
    %     for i=1:n_bands1
    %         if f1(i,1)<f1(i,2) % First filter definition mode [lambda_min, lambda_max, transmission]
    %             w_ch1 = find((l >= f1(i,1) & l <= f1(i,2)));
    %         else % Second filter definition mode [lambda_central, bandpass, transmission]
    %             w_ch1 = find((l >= (f1(i,1)-f1(i,2)/2)) & (l <= (f1(i,1)+f1(i,2)/2)));
    %         end
    %         if ~isempty(w_ch1)
    %             n1 = n1+fix(sum(s(w_ch1))*f1(i,3));
    %         end
    %     end
    n1=sum(s.*filter_profiles.ch1);

    %Apply Poissonian law if needed
    if im_par.stochastic_spectrum_off==1
        n1=poissrnd(n1);
    end

    %Second channel
    %     for i=1:n_bands2
    %         if f1(i,1)<f1(i,2) % First filter definition mode [lambda_min, lambda_max, transmission]
    %             w_ch2 = find((l >= f2(i,1) & l <= f2(i,2)));
    %         else
    %             w_ch2 = find((l >= (f2(i,1)-f2(i,2)/2)) & (l <= (f2(i,1)+f2(i,2)/2)));
    %         end % Second filter definition mode [lambda_central, bandpass, transmission]
    %         if ~isempty(w_ch2)
    %             n2 = n2+fix(sum(s2(w_ch2))*f2(i,3));
    %         end
    %     end
    n2=sum(s2.*filter_profiles.ch2);

    %Apply Poissonian law if needed
    if im_par.stochastic_spectrum_off==1
        n2=poissrnd(n2);
    end

end

if im_par.debug
    disp(['# of photons emitted in channel 1: ',num2str(n1)]);
    if im_par.two_channel==1
        disp(['# of photons emitted in channel 2: ',num2str(n2)]);
    end
end

end

