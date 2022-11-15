function [x_ch,y_ch,lambda_mean] = get_filter_profile(ch_filter)

% NAME:
%	get_filter_profile
%
% PURPOSE:
%       Read the profile of a fluorescence filter
%
% INPUTS:
%   ch_filter: parameters for the filter ([lambda_min, lambda_max,
%   transmission]), possibly multiband
%
% OUTPUTS:
%	x_ch
%   y_ch
%   lambda_mean: the mean wavelength over the filter bandpass
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2017.
%	D.Bourgeois, September 2022. Add filter definition [lambda_min, lambda_max, transmission]
%-

max_lambda=1000; % We limit to 1000 nm the microscope bandpss

if size(ch_filter,2) ~= 3
    warning('Filter definition must have 3 parameters !');
    return
end
n_bands=size(ch_filter,1); % # of bands in the filter

%check filter
xs_ch=zeros(1,n_bands);
xe_ch=zeros(1,n_bands);
lm=zeros(1,n_bands);
for i=1:n_bands
    if ch_filter(i,1)<ch_filter(i,2) % Format is lambda_min, lambda_max, transmission
        xs_ch(i)=fix(ch_filter(i,1));
        xe_ch(i)=ceil(ch_filter(i,2));
    else % Format is lambda_central, bandpass, transmission
        xs_ch(i)=fix(ch_filter(i,1)-ch_filter(i,2)/2);
        xe_ch(i)=fix(ch_filter(i,1)+ch_filter(i,2)/2);
    end
    lm(i)=0.5*(xs_ch(i)+min([xe_ch(i),max_lambda]));
end

lambda_mean=round(mean(lm));

if n_bands>1
    for i=2:n_bands
        if xs_ch(i)<xe_ch(i-1)
            warning('Filter bands must be in increasing order !');
            return
        end
    end
end

x_ch=[];
y_ch=[];
for i=1:n_bands
    x_ch_loc=[xs_ch(i)-1, xs_ch(i),xe_ch(i), xe_ch(i)+1]; % to get the filter edges
    y_ch_loc=[0, 1, 1, 0]*ch_filter(i,3);
    
    if isinf(xe_ch(i))
        x_ch_loc=x_ch_loc(1:3);
        x_ch_loc(end)=max_lambda;
        y_ch_loc=y_ch_loc(1:3);
    end

    if ch_filter(i,2)==0
        y_ch_loc(:)=0;
    end
    x_ch=horzcat(x_ch,x_ch_loc);
    y_ch=horzcat(y_ch,y_ch_loc);
end
    
    
