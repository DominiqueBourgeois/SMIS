function y = get_exp_aniso_model(a, t)
% PURPOSE:
%
%       For fitting, returns a (stretched) exponential decay model taking
%       into account anisotropy of single-molecules orientations
%
% CATEGORY:
%       Data processing
%
% CALLING SEQUENCE:
%       y = get_exp_aniso_model(a, t)
%
% INPUTS:
%       a: the variables
%       t: the time data
%
%
% MODIFICATION HISTORY:
%       D.Bourgeois, November 2012.

n=1000; % Create distribution of theta angles
theta=pi*(0:n)/n-pi/2;
p=cos(theta); 
p_norm=p/sum(p); % Probability of a randomly oriented molecule to have this theta
rate_change=3/2*p.^2; 

y=zeros(length(t),1);
for k=1:length(t)
    y(k)=a(1)*sum(p_norm.*rate_change.*exp(-rate_change*(a(2)*a(4)*t(k)).^a(5)))+a(3);
end
y=y*a(1)/y(1);



