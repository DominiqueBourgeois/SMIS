function y = get_exp_aniso_linear_model(a, t)
% PURPOSE:
%
%       For fitting, returns a (stretched) exponential decay model taking
%       into account anisotropy of single-molecules orientations
%
% CATEGORY:
%       Data processing
%
% CALLING SEQUENCE:
%       y = get_exp_aniso_linear_model(a, t)
%
% INPUTS:
%       a: the variables
%       t: the time data
%
%
% MODIFICATION HISTORY:
%       D.Bourgeois, November 2012.

n=64; % Create distribution of theta angles
theta=pi*(0:n)/n-pi/2;

m=64; % Create distribution of phi angles
phi=2*pi*(0:m)/m;

p_theta_phi=zeros(length(theta),length(phi)); % Probability of a randomly oriented molecule to have certain phi & theta
k_theta_phi=zeros(length(theta),length(phi)); % Excitation rates at phi & theta


for k=1:length(theta)
    p_theta_phi(k,:)= cos(theta(k));
    k_theta_phi(k,:)= 3*cos(theta(k))^2*cos(phi).^2; % Factor of 3 = 3/2 (for theta) * 2 (for phi)
end

p_theta_phi=p_theta_phi/sum(sum(p_theta_phi)); %Normalize to 1

y=zeros(length(t),1);

for k=1:length(t)
    y(k)=a(1)*sum(sum(p_theta_phi.*k_theta_phi.*exp(-k_theta_phi*(a(2)*a(4)*t(k)).^a(5))));
end
y=y*a(1)/y(1)+a(3);



