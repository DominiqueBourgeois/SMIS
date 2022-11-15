function [s1, s2] = filter_em_spectrum_ensemble(spectrum, filters, im_par)
% NAME:
%	filter_em_spectrum_ensemble
%
% PURPOSE:
%	Filter spectrum by 2 filters for 2-channel ensemble experiments
%
% CATEGORY:
%	Signal, image processing.
%
% CALLING SEQUENCE:
%	[fr1, fr2] = filter_em_spectrum_ensemble(spectrum, par)
%
% INPUTS:
%   spectrum: the spectrum to be filtered
%   filters: the channel filter profiles
%	im_par: parameters containing filters infos
%
% OUTPUTS:
%	s1, s2:  signal emitted through each filter
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
%	D.Bourgeois, September 2021.


s=spectrum(:,2);

%Channel 1
f=filters.ch1; % Filter profile for ch1

if im_par.two_channel==1 && im_par.filters.add_dichroic==1
    df=filters.df; % The dichroic
    fs=s.*f.*(1-df);
else
    fs=s.*f;
end

s1=sum(fs);

%Channel 2
if im_par.two_channel==1
    f=filters.ch2; % Filter profile for ch2

    if im_par.filters.add_dichroic==1
        df=filters.df; % The dichroic
        fs=s.*f.*df;
    else
        fs=s.*f;
    end
    s2=sum(fs);
else
    s2=0;
end



% f1 = im_par.filters.ch1_filter; % Get filter 1
% if im_par.two_channel==1; f2 = im_par.filters.ch2_filter; end % Get filter 2

% n_bands1=size(f1,1); % # of bands in the filter
% if im_par.two_channel==1; n_bands2=size(f2,1); end % # of bands in the filter

% s1=0;
% for i=1:n_bands1
%     w_ch1 = find((spectrum(:,1) >= (f1(i,1)-f1(i,2)/2)) & (spectrum(:,1) <= (f1(i,1)+f1(i,2)/2)));
%     if ~isempty(w_ch1)
%         s1 = s1+sum(spectrum(w_ch1,2))*f1(i,3);
%     end
% end
%
% s2=0;
% if im_par.two_channel==1
%     for i=1:n_bands2
%         w_ch2 = find((spectrum(:,1) >= (f2(i,1)-f2(i,2)/2)) & (spectrum(:,1) <= (f2(i,1)+f2(i,2)/2)));
%         if ~isempty(w_ch2)
%             s2 = s2+sum(spectrum(w_ch2,2))*f2(i,3);
%         end
%     end
% end


