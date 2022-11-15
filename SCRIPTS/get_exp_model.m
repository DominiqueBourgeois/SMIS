function y = get_exp_model(a, t)
%+
% NAME:
%       GET_MODEL
%
% PURPOSE:
%
%       For fitting, returns a stretched exponential decay model
%
% CATEGORY:
%       Data processing
%
% CALLING SEQUENCE:
%       y = get_exp_model(a, t)
%
% INPUTS:
%       a: the variables
%       t: the time data
%
%
% MODIFICATION HISTORY:
%       D.Bourgeois, November 2012.


y=a(1)*exp(-(a(2)*a(4)*t).^a(5))+a(3);




