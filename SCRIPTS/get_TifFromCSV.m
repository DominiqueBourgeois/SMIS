function im=get_TifFromCSV(filename, par)

% NAME:
%	get_TifFromCSV
%
% PURPOSE:
%       Read a .CSV file containing coordinates and create a .TIF image or
%       stack. The columns names should contains 'x' or 'X', 'y' or 'Y',
%       'z' or 'Z' (for 3D samples), and possibly 'ids' (for pattern ids)
% INPUTS:
%   par: image parameters
%
% OUTPUTS:
%	im : the image or kernel
%
% MODIFICATION HISTORY:
%	D.Bourgeois, May 2015.
%	D.Bourgeois, March 2023: added read tif files for 3D images.
%	D.Bourgeois, June 2025: Add possibility to read .csv virtual samples containing sample coordinates [nm]
%-


% Read the CSV file
data = readtable(filename);

% define image border
image_border = 10;

% Check that the right column names exist
read_ok=1;
columnNames = data.Properties.VariableNames;
hasX = any(strcmpi(columnNames, 'x'));
hasY = any(strcmpi(columnNames, 'y'));
hasZ = any(strcmpi(columnNames, 'z'));
hasids = any(strcmpi(columnNames, 'ids'));

R = par.raster/par.binning; % Resolution of the produced .tif image

% Display a warning if either 'x' or 'y' or 'z' columns are not found
if par.simul_3D==1
    if ~hasX || ~hasY || ~hasZ
        warningMessage = 'Warning: The CSV file does not contain the required columns.';
        if ~hasX
            warningMessage = [warningMessage ' Column ''x'' or ''X'' is missing.'];
        end
        if ~hasY
            warningMessage = [warningMessage ' Column ''y'' or ''Y'' is missing.'];
        end
        if ~hasZ
            warningMessage = [warningMessage ' Column ''z'' or ''Z'' is missing.'];
        end
        read_ok=0;
    end
else
     if ~hasX || ~hasY
        warningMessage = 'Warning: The CSV file does not contain the required columns.';
        if ~hasX
            warningMessage = [warningMessage ' Column ''x'' or ''X'' is missing.'];
        end
        if ~hasY
            warningMessage = [warningMessage ' Column ''y'' or ''Y'' is missing.'];
        end
        read_ok=0; 
     end
end
if read_ok==0
    im=[];
    warningMessage = [warningMessage ' Sample cannot be loaded.'];
    MyDlg=warndlg(warningMessage, 'Column Check Warning');
    waitfor(MyDlg)
    return
end


% Extract X, Y, Z columns
if ismember('x', data.Properties.VariableNames)
    X = data.x; 
elseif ismember('X', data.Properties.VariableNames)
    X = data.X; 
end
if ismember('y', data.Properties.VariableNames)
    Y = data.y; 
elseif ismember('Y', data.Properties.VariableNames)
    Y = data.Y; 
end
if par.simul_3D==1
    if ismember('z', data.Properties.VariableNames)
        Z = data.z; 
    elseif ismember('Z', data.Properties.VariableNames)
        Z = data.Z; 
    end
end

% Convert coordinates to indices based on resolution R
X_idx = round(X / R);
Y_idx = round(Y / R);
if par.simul_3D==1
    Z_idx = round(Z / R);
end

if hasids
    ids = data.ids; % Pattern ids
    disp(['Number of subpattern ids found = ',num2str(numel(unique(ids)))]);
else
    ids = 1;
    disp('No subpattern id found, setting all coordinates to id = 1');
end


% Make sure coodinates are positive
X_idx=X_idx-min(X_idx)+1+image_border;
Y_idx=Y_idx-min(Y_idx)+1+image_border;
if par.simul_3D==1
    Z_idx=Z_idx-min(Z_idx)+1+image_border;
end

% Determine the size of the 3D matrix and initialize matrix

maxX = max(X_idx);
maxY = max(Y_idx);
if par.simul_3D==1
    maxZ = max(Z_idx);
    im = zeros(maxY+image_border, maxX+image_border, maxZ+image_border);
    im(sub2ind(size(im), Y_idx, X_idx, Z_idx)) = ids;
else
    im = zeros(maxY+image_border, maxX+image_border);
    % Set the voxels corresponding to the points to ids values
    im(sub2ind(size(im), Y_idx, X_idx)) = ids;
end

