function psf_em = get_sm_emitted_psf(n_phot,pos,psf, par)

% NAME:
%	get_sm_emitted_psf
%
% PURPOSE:
%       Simulates a 2D image created by a single fluorophore
%       through a microscope with a given PSF. If N photons are emitted,
%       these N photons end up in the various pixels of the normalized PSF
%       image so that p(photon ends up in pixel i)= PSF(i)
%
% CATEGORY:
%	Super-resolution
%
% CALLING SEQUENCE:
%	s_em = get_sm_emitted_psf(n_phot,pos,psf, im_par)
%
% INPUTS:
%	n_phot: the number of photons in the peak
%   pos: position of the molecule on the image
%   psf: the psf image
%   im_par: general paramters 
%
% OUTPUTS:
%	peak : The 2D peak
%
% MODIFICATION HISTORY:
%	D.Bourgeois, May 2015. December 2017
%-

% Shift, normalize and multiply psf by 1000
psf_s = size(psf);
scale=1000;

% [xi,yi] = meshgrid((1:psf_s(1))+pos(1),(1:psf_s(2))+pos(2));
% [xi,yi] = meshgrid((1:psf_s(1))-pos(1),(1:psf_s(2))-pos(2)); % Corrected bug 20/12/2017 shift should be in negative direction
[xi,yi] = meshgrid((1:psf_s(2))-pos(1),(1:psf_s(1))-pos(2)); % Corrected bug 20/12/2017 shift should be in negative direction

% psf_i=interp2(psf,xi,yi,'spline'); % Could be faster with linear interpolation
psf_i=interp2(psf,xi,yi,'linear',0); % Could be faster with linear interpolation

psf_n=scale*psf_i/sum(psf_i(:));

%indices of emitted phot
i_em = scale*rand(n_phot,1);

%define sm psf 
psf_em=zeros(psf_s);

%get the cumulative integral of the input psf
psf_c=zeros(psf_s);
psf_c(1)=psf_n(1);
for i=2:length(psf_c(:))
    psf_c(i)=psf_c(i-1)+psf_n(i);
end

%look where the i_em end up on psf_c: i_em(i) must be associated to the
%smallest psf_c value which is > i_em(i).
for i=1:length(i_em)
    sel_i = find(psf_c >= i_em(i),1,'first');
    psf_em(sel_i)=psf_em(sel_i)+1;
end

psf_em=transpose(psf_em); % Added 12/2017 was inconsistent before

if par.debug
    subplot(2,1,1);
    plot(psf_n);
    subplot(2,1,2);
    plot(psf_em);
end

