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
%	D.Bourgeois, March 2023: added read tif files for 3D images.
%-

[~,~,ext]=fileparts(filename);

if par.simul_3D==0 && (strcmp(ext,'.tif')==1 || strcmp(ext,'.tiff')==1)
    a=imread(filename);
elseif par.simul_3D==1 && strcmp(ext,'.mat')==1
    load(filename,'kernel_pattern');
    a=int16(kernel_pattern);
elseif par.simul_3D==1 && (strcmp(ext,'.tif')==1 || strcmp(ext,'.tiff')==1)
    [kernel_pattern, ~] = read_SMIS_MultipageTiff(filename);
    a=int16(kernel_pattern);
else
    disp('Error: Unaccepted pattern image format: use Tif or Matlab format !');
    return;
end
if par.simul_3D==1
    disp(['Number of slices in 3D image: ',num2str(size(a,3))]);
end
[a,n,m,z]=check_im_size(a, par.binning);

