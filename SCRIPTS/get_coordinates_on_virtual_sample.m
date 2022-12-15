function [x,y,z]=get_coordinates_on_virtual_sample(x0,y0,z0, im_par)

%
% PURPOSE:
%   Get the coordinates of single molecules on the virtual sample detector from the
%   coordinates on the detector (x0, y0).
%
% INPUTS:
%   x0, y0, z0: the single molecules (with coordinates on detector) in raster units
%	im_par: the imaging parameters
%
% OUTPUTS:
%   x, y, z: the coordinates in raster units on the virtual sample
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2022

% For an image that goes from 1 to N in indices, x and y go
% from 0.5 to N+0.5. The xy value at center of pixel (i,j) is
% x=i and y=j. The xy value throughout pixel (i,j) go from
% x=i-0.5; y=j-0.5 to x=i+0.5; y=j+0.5.
x=(x0-0.5)*im_par.binning+0.5;
y=(y0-0.5)*im_par.binning+0.5;

if im_par.simul_3D==1
    z=(z0-0.5)*im_par.binning+0.5;
else
    z=[];
end

