function R0 = get_R0(exc, em, eps, QY, par, kappa2)
% NAME:
%   get_sm_emitted_spectrum
%
% PURPOSE:
%	%This program calculate Forster distance (R0) using Lakowicz eqn(13.5)
%
% CATEGORY:
%	Signal, image processing.
%
% CALLING SEQUENCE:
%	r0 = get_R0(sm_par,sm_par2, par, kappa)
%
% INPUTS:
%   sm_par: general parameters for the donor
%   sm_par_2: general parameters for the acceptor
%   par: other general parameters
%   kappa: orientation factor
%
% OUTPUTS:
%	r0: the Forster distance
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
%   MSG 11-24-2009
%	D.Bourgeois, April 2011.
%

avogadro = 6.02e+23; % [mol-1]

% exc=sm_par_2.exc_spectrum; %the excitation spectrum of the acceptor
% em=sm_par.em_spectrum; %the emission spectrum of the donor
% eps=sm_par_2.extinction_coeff;
% QY=sm_par.quantum_yield;  %donor Quantum Yield
if nargin < 6
    kappa2=2/3; %Kappa^2 assumed to be 2/3 for random distribtuion
end

n=1.4;%refractive index for biomolecules;

Fd=em(:,2); %Normalized fluorescence emission of donor
wl_em=em(:,1);

Ea=eps*exc(:,2); %Fluorescence excitation of acceptor
wl_exc=exc(:,1);

%the two spectra must have the same number of elements
Ea2=Fd;
for i=1:length(wl_em)
    if i==1
        step=wl_em(i+1)-wl_em(i);
    else
        step=wl_em(i)-wl_em(i-1);
    end
    tmp_var=abs(wl_exc-wl_em(i));
    w=find(tmp_var==min(tmp_var)); %find nearest point in excitation spectrum   
    w=w(1);
    if min(tmp_var)<step % an element of acceptor is close enough
        Ea2(i)=Ea(w);
    else
        Ea2(i)=0;
    end
end

%Do integral 
I=sum(Fd.*Ea2.*wl_em.^4)/sum(Fd); 

%epsilon in Ea2 is in L.mol-1.cm-1 => dm^3.cm-1 => 10^27 A^3/(10^8 A) => 10^19 A^2
%lambda in nm => 10^4A
% so 10^23 A
%Then in Lakowski the real value is 9*log(10)/(128*pi^5), not 9000*log(10)/(128*pi^5) 
R0=(1e+23*9*log(10)/(128*pi^5)/avogadro)^(1/6)*(kappa2*n^(-4)*QY*I)^(1/6);
% This agrees with Siyath expression R0=0.211*(kappa2*n^(-4)*qd*I)^(1/6); %R0 in Å

if par.debug
    plot(wl_em,Fd*eps,'r');
    hold on
    plot(wl_exc,Ea,'b');
    plot(wl_em,Ea2,'b');
    hold off
end
end