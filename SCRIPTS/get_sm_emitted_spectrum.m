function s_em = get_sm_emitted_spectrum(n_phot,em_spectrum,im_par)
% NAME:
%   get_sm_emitted_spectrum
%
% PURPOSE:
%	Get emitted spectrum by a single molecule
%
% CATEGORY:
%	Signal, image processing.
%
% CALLING SEQUENCE:
%	s_em = get_sm_emitted_spectrum(n_phot,em_spectrum, im_par);
%
% INPUTS:
%   n_phot: # of emitted photons
%	em_spectrum: the reference emission in solution.
%   im_par: general parameters
%
% OUTPUTS:
%	s_em: the emitted spectrum by the single molecule
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
%	D.Bourgeois, April 2011.
%	D.Bourgeois, May 2015. Corrected small bug.

scale=1000;
s=em_spectrum(:,2);
s=s*scale/sum(s); %normalize the spectrum to 'scale'.
ss = size(s);

%Apply Poissonian law
n_phot=poissrnd(n_phot);

%indices of emitted phot
i_em = scale*rand(n_phot,1);

%define sm emmission spectrum
s_em=zeros(ss(1),2);
s_em(:,1)=em_spectrum(:,1);

%get the cumulative integrale
sc=zeros(ss);
sc(1)=s(1);
for i=2:length(sc)
    sc(i)=sc(i-1)+s(i);
end

%look where the i_em end up on sc
for i=1:length(i_em)
    sel_i = find(sc >= i_em(i),1,'first');
    s_em(sel_i,2)=s_em(sel_i,2)+1;
end

if im_par.debug
    plot(em_spectrum(:,1),em_spectrum(:,2)*max(s_em(:,2))/max(em_spectrum(:,2)));
    hold on
    bar(s_em(:,1),s_em(:,2));
    hold off
end

