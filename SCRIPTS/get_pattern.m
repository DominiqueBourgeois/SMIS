function [a,n,m,z] = get_pattern(filename, par)

% NAME:
%	get_pattern
%
% PURPOSE:
%       Read a pattern image (.tif) or 3D-kernal (.mat)
% INPUTS:
%   par: image parameters
%
% OUTPUTS:
%	a : the image or kernel
%   n,m,z : the dimensions [raster]
%
% MODIFICATION HISTORY:
%	D.Bourgeois, May 2015.
%-

[~,~,ext]=fileparts(filename);
if strcmp(ext,'.tif')==1 || strcmp(ext,'.tiff')==1
    a=imread(filename);
elseif strcmp(ext,'.mat')==1
    load(filename,'kernel_pattern');
    a=int16(kernel_pattern);
else
    disp('Error: Unaccepted pattern image format: use Tif or Matlab format !');
    return;
end
[a,n,m,z]=check_im_size(a, par.binning);
