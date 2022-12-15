function [x_df,y_df,lambda_mean] = get_dichroic_filter_profile(df)

% NAME:
%	get_dichroic_filter_profile
%
% PURPOSE:
%       Read the profile of a dichroic filter
%
% INPUTS:
%   df: parameters for the filter ([lambda_cutoff, bandpass,
%   transmission_lowpass, transmission_highpass])
%
% OUTPUTS:
%	x_df
%   y_df
%   lambda_mean: the mean wavelength over the filter bandpass
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2022
%-

if size(df,2) ~= 4
    warning('Filter definition must have 4 parameters !');
    return
end

%Plot limits
min_l=200;
max_l=800; 

if min_l<(df(1)-df(2)/2) && max_l>(df(1)+df(2)/2) % General case
    x_df=[min_l,df(1)-df(2)/2,df(1)+df(2)/2, max_l];
    y_df=[df(3), df(3), df(4),df(4)];
elseif min_l<(df(1)-df(2)/2) && max_l<=(df(1)+df(2)/2)
    x_df=[min_l,df(1)-df(2)/2,df(1)+df(2)/2];
    y_df=[df(3), df(3), df(4)];
elseif min_l>=(df(1)-df(2)/2) && max_l>(df(1)+df(2)/2)
    x_df=[df(1)-df(2)/2,df(1)+df(2)/2,max_l];
    y_df=[df(3), df(4), df(4)];
else
    x_df=[df(1)-df(2)/2,df(1)+df(2)/2];
    y_df=[df(3), df(4)];
end

lambda_mean=round(df(1));

