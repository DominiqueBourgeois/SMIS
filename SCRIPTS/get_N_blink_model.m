function y = get_N_blink_model(a, x)
%+
% NAME:
%       GET_N_BLINK_MODEL
%
% PURPOSE:
%
%       For fitting, returns geometric law: y = eta^x*(1-eta)
%
%       D.Bourgeois, July 2013.

a0=a(1);
eta=a(2);

y=a0*(eta.^x)*(1-eta);




