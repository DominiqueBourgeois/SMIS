function [x_ch,y_ch,lambda_mean] = get_exp_filter_profile(filename)

% NAME:
%	get_exp_filter_profile
%
% PURPOSE:
%       Read the profile of a fluorescence filter from commercial filter
%
% INPUTS:
%   filename: text file containing x and y columns with wavelength and
%   transmission values
%
% OUTPUTS:
%	x_ch
%   y_ch
%   lambda_mean: the mean wavelength over the filter bandpass
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2022.
%-

% We limit to 200-1000 nm the microscope bandpss
lambda_min=200;
lambda_max=1000; 

try
    data=table2array(readtable(filename));
    x=data(:,1);
    y=data(:,2);
    %Cut above 1000 nm and below 200 nm
    y=y(x>=lambda_min & x<=lambda_max);
    x=x(x>=lambda_min & x<=lambda_max);
    %Get data every nm
    min_lambda=fix(min(x));
    max_lambda=ceil(max(x));
    x_ch=min_lambda:max_lambda;
    y_ch=interp1(x,y,x_ch,"linear");
    lambda_mean=sum(y_ch.*x_ch)/sum(y_ch);
catch
    warndlg('Wrong file format !','Warning')
    x_ch=[];
    y_ch=[];
    lambda_mean=[];
end
