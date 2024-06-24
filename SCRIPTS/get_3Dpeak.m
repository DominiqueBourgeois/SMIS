function peak=get_3Dpeak(x,y,z,n_phot,sm_par,im_par,state,channel)
% NAME:
%	GET_3DPEAK
%
% PURPOSE:
%       Simulates a z-dependent 2D image created by a single fluorophore
%       through a microscope with a z-dependent PSF and shot-noise
% INPUTS:
%   x, y, z: coordinates of the single molecule on the detector image [raster
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
%	D.Bourgeois, September 2020: Returns empty peak if outside z range. Changed 3D case
%	D.Bourgeois, March 2022: Correct mapping between z position and psf slice
%-

show_peak=0; % set to 1 to see the peak

if channel==1
    psf=sm_par.psf_par_ch1(state).psf;
else
    psf=sm_par.psf_par_ch2(state).psf;
end

% take the factional xy coordinates
% For an image that goes from 1 to N in indices, x and y go
% from 0.5 to N+0.5. The xy value at center of pixel (i,j) is
% x=i and y=j. The xy value throughout pixel (i,j) go from
% x=i-0.5; y=j-0.5 to x=i+0.5; y=j+0.5.
r_xy = [x-round(x),y-round(y)]; % rounded x and y

%get_3Dpeak
%Get the psf slice adequate for z position

% z goes from 0.5 to im_par.nz+0.5
% n_psf_slice must go from 1 to im_par.psf_n_zslices, ie f(0.5)=1 & f(im_par.nz+0.5)=im_par.psf_n_zslices
a_=(im_par.psf_n_zslices-1)/(im_par.nz);
b_ = 1-0.5*a_;

% n_psf_slice=(z+0.5)*im_par.psf_n_zslices/im_par.nz;
n_psf_slice=z*a_+b_;

bottom_psf_slice=max([fix(n_psf_slice),1]);
top_psf_slice=min([ceil(n_psf_slice),im_par.psf_n_zslices]);
if bottom_psf_slice<=size(psf,2) && top_psf_slice<=size(psf,2) && bottom_psf_slice>0 && top_psf_slice>0 % Check that molecule has not diffused away
    if mod(n_psf_slice,bottom_psf_slice)>0  %get some kind of interpolation !
        % Handle case where psf at top_psf_slice and bottom_psf_slice are
        % not of same size
        if psf(top_psf_slice).rbox_x==psf(bottom_psf_slice).rbox_x && psf(top_psf_slice).rbox_y==psf(bottom_psf_slice).rbox_y
            psf_slice=(n_psf_slice-bottom_psf_slice)*psf(top_psf_slice).plane+ (top_psf_slice-n_psf_slice)*psf(bottom_psf_slice).plane;           
        else
            rx=max([psf(bottom_psf_slice).rbox_x,psf(top_psf_slice).rbox_x]);
            ry=max([psf(bottom_psf_slice).rbox_y,psf(top_psf_slice).rbox_y]);
            loc_psf_top=zeros(2*ry+1,2*rx+1);
            loc_psf_bottom=zeros(2*ry+1,2*rx+1);
            loc_psf_top(ry-psf(top_psf_slice).rbox_y+1:ry+psf(top_psf_slice).rbox_y+1,rx-psf(top_psf_slice).rbox_x+1:rx+psf(top_psf_slice).rbox_x+1)=psf(top_psf_slice).plane;
            loc_psf_bottom(ry-psf(bottom_psf_slice).rbox_y+1:ry+psf(bottom_psf_slice).rbox_y+1,rx-psf(bottom_psf_slice).rbox_x+1:rx+psf(bottom_psf_slice).rbox_x+1)=psf(bottom_psf_slice).plane;
            psf_slice=(n_psf_slice-bottom_psf_slice)*loc_psf_top + (top_psf_slice-n_psf_slice)*loc_psf_bottom;           
        end
    else
        psf_slice=psf(bottom_psf_slice).plane;
    end
    peak=get_sm_emitted_psf(n_phot,r_xy, psf_slice, im_par);
else % if molecule out of z range, do not show it
    peak=[];
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
end

