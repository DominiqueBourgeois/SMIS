function peak=get_2Dpeak(x,y,n_phot,sm_par,im_par,state,channel)

% NAME:
%	GET_2DPEAK
%
% PURPOSE:
%       Simulates a 2D image created by a single fluorophore
%       through a microscope with a 2D Gaussian PSF and shot-noise
%       Note that we do not use the .psf defined in get_psf.m. This could
%       be done like in 3D case using get_sm_emitted_psf.m, to be checked
%       if this is computationally faster or not.
% INPUTS:
%   x, y: coordinates of the single molecule on the detector image [raster
%   units]
%	n_phot: the number of photons in the peak
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%   state: number of the fluorescence state
%   channel: detector channel number
%
% OUTPUTS:
%	peak : The 2D peak
%
% MODIFICATION HISTORY:
%	D.Bourgeois, April 2011. September 2019
%-

show_peak=0; % set to 1 to see the peak

if channel==1
    fwhm=sm_par.psf_par_ch1(state).psf_width;
    rbox=sm_par.psf_par_ch1(state).rbox;
else
    fwhm=sm_par.psf_par_ch2(state).psf_width;
    rbox=sm_par.psf_par_ch2(state).rbox;
end

%fwhm of the psf in pixels
fwhm_r = fwhm/im_par.raster;

%sigma of the psf
sig=fwhm_r/2.35;

% take the factional xy coordinates
% For an image that goes from 1 to N in indices, x and y go
% from 0.5 to N+0.5. The xy value at center of pixel (i,j) is
% x=i and y=j. The xy value throughout pixel (i,j) go from
% x=i-0.5; y=j-0.5 to x=i+0.5; y=j+0.5.
r_xy = [x-round(x),y-round(y)]; % rounded x and y
pos=[rbox+1,rbox+1] + r_xy ; % position in box

%distribute these n photons throughout the PSF
sx=2*rbox+1;
sy=2*rbox+1;
peak=zeros(sx,sy);
rn=sig*randn(n_phot,2); % 2D Gaussian distribution of width = sig
xs=round(rn(:,1)+pos(1)); %photons will end up in most adjacent pixels
ys=round(rn(:,2)+pos(2));
for k=1:n_phot
    x=xs(k);
    y=ys(k);
    if (x > 0) && (x <= sx) && (y > 0) && (y <= sy)
        peak(x,y) = peak(x,y)+1;
    end
end

if show_peak==1
    figure(1)
    surf(peak);
    xlabel('X [raster]','fontsize',8,'fontweight','b')
    ylabel('Y [raster]','fontsize',8,'fontweight','b')
    title(gca,'Single Molecule Peak','FontWeight','bold');
    disp(['# of emitted photons in each fluorescent state:' , num2str(n_phot)]);
    input('Ok ? ','s');
end


