function [x2,y2,z2, sm] = get_linkage_error_pct(x,y,z,sm, im_par, sm_par)
%
% NAME:
%	get_linkage_error
% PURPOSE:
%
% INPUTS:
%   x, y, z: coordinate of the single molecule on detector image [raster units]
%   sm: the single molecule (with coordinates on high-resolution image) in raster units
%	im_par: the imaging parameters
%	sm_par: the sm parameters
%
% OUTPUTS:
%	x2,y2,z2: the corrected x y z for linkage error [raster units]
%   sm: The updated single molecule for linkage error
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2022

% Get the proper indices
lx_idx=1;
ly_idx=2;
lz_idx=3;
le_set_idx=4;

if sm_par.linkage_fixed==1 % Case of fixed linkage error through data collection
    if isempty(sm{lx_idx}) % Set linkage error once
        theta=2*pi*rand; % Pick random angle
        l=(sm_par.linkage_length+randn*sm_par.linkage_std)/im_par.raster; % Length of linkage error [raster unit]
        if im_par.simul_3D==1
            phi=pi*rand; % Pick random angle
            lx=cos(theta)*sin(phi)*l;
            ly=sin(theta)*sin(phi)*l;
            lz=cos(phi)*l;
        else
            lx=cos(theta)*l;
            ly=sin(theta)*l;
            lz=[];
        end
    else % Just use the fixed linkage error
        lx=sm{lx_idx};
        ly=sm{lx_idx};
        lz=sm{lx_idx};
    end
else % Case of random linkage error through data collection
    if sm{le_set_idx}==im_par.current_frame % The linkage error has been already calculated for the other channel
        lx=sm{lx_idx};
        ly=sm{lx_idx};
        lz=sm{lx_idx};
    else
        theta=2*pi*rand; % Pick random angle
        l=(sm_par.linkage_length+randn*sm_par.linkage_std)/im_par.raster; % Length of linkage error
        if im_par.simul_3D==1
            phi=pi*rand; % Pick random angle
            lx=cos(theta)*sin(phi)*l;
            ly=sin(theta)*sin(phi)*l;
            lz=cos(phi)*l;
        else
            lx=cos(theta)*l;
            ly=sin(theta)*l;
            lz=[];
        end
    end
end

% Keep the linkage errors
sm{lx_idx}=lx;
sm{ly_idx}=ly;
sm{lz_idx}=lz;
sm{le_set_idx}=im_par.current_frame;

% Assign the linkage errors
x2=x+lx;
y2=y+ly;
if im_par.simul_3D==1
    z2=z+lz;
else
    z2=[];
end



