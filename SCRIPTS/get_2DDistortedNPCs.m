function im=get_2DDistortedNPCs(par)

% NAME:
%	get_2DDistortedNPCs
%
% PURPOSE:
%       Create a 2D image from NPC coordinates
% INPUTS:
%   par: image parameters
%
% OUTPUTS:
%	im : the image or kernel
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2025
%-


% Read the CSV file
data = par.data;

R = par.raster/par.binning; % Resolution of the produced .tif image

% Display a warning if either 'x' or 'y' or 'z' columns are not found

% Extract X, Y, Z columns
X = data.x;
Y = data.y;
% Convert coordinates to indices based on resolution R
X_idx = round(X / R);
Y_idx = round(Y / R);

if hasids
    ids = data.ids; % Pattern ids
    disp(['Number of subpattern ids found = ',num2str(numel(unique(ids)))]);
else
    ids = 1;
end

% Make sure coodinates are positive
X_idx=X_idx-min(X_idx)+1;
Y_idx=Y_idx-min(Y_idx)+1;

% Determine the size of the 3D matrix and initialize matrix

maxX = max(X_idx);
maxY = max(Y_idx);
if simul_3D==1
    maxZ = max(Z_idx);
    im = zeros(maxX, maxY, maxZ);
    im(sub2ind(size(im), X_idx, Y_idx, Z_idx)) = ids;
else
    im = zeros(maxX, maxY);
    % Set the voxels corresponding to the points to ids values
    im(sub2ind(size(im), X_idx, Y_idx)) = ids;
end

