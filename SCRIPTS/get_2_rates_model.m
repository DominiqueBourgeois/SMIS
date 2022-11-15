function y = get_2_rates_model(a, t)
%+
% NAME:
%       get_2_rates_model
%
% CALLING SEQUENCE:
%       y = get_2_rates_model(a, t)
%
% INPUTS:
%       a: the variables
%       t: the time data
%
%
% MODIFICATION HISTORY:
%       D.Bourgeois, Janvier 2015.

a0_1=a(1);
a0_2=a(2);
offset=a(3);
k1=a(4);
k2=a(5);

y=a0_1*exp(-k1*t)+a0_2*exp(-k2*t)+offset;




