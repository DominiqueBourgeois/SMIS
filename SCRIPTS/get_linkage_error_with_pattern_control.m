function [x2,y2,z2, sm, ok] = get_linkage_error_with_pattern_control(x,y,z,sm, im_par, sm_par)
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
%   Not implemented currently: to ensure that the linkage distance between
%   actual target and fluorophore induces a change of subpattern (eg target
%   on a membrane and fluorophore sticking out in external medium), the
%   resolution of the virtual sample image should be on the nm range, ie
%   binning factors of ~100. Coding not finished for 3D case.

% Get the proper indices
c_sp_idx=1;
lx_idx=2;
ly_idx=3;
lz_idx=4;
ok=1; % Default value

if sm_par.linkage_pattern_sp_id==0 % Case of linkage error forcing the fluorophore to stick out of current pattern
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
            % Keep these values
            sm{lx_idx}=lx;
            sm{ly_idx}=ly;
            sm{lz_idx}=lz;
        else % Just use the fixed linkage error
            lx=sm{lx_idx};
            ly=sm{lx_idx};
            lz=sm{lx_idx};
        end
    else % Case of random linkage error through data collection
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
else
    sz=im_par.binning*[im_par.n, im_par.m]; % Size of virtual sample image
    max_trials=1000; % Maximun number of trials
    if sm_par.linkage_fixed==1 % Case of fixed linkage error through data collection
        if isempty(sm{lx_idx}) % Set linkage error once
            l=(sm_par.linkage_length+randn*sm_par.linkage_std)/im_par.raster; % Length of linkage error [raster unit]
            stick_out=0;
            n_trial=1;
            w=sm_par.w_patterns(sm_par.n_sp_id==sm{c_sp_idx}).w; % Indices on virtual sample of the sm current subpattern
            while stick_out==0 && n_trial<max_trials
                theta=2*pi*rand; % Pick random angle
                if im_par.simul_3D==1
                    phi=pi*rand; % Pick random angle
                    lx=cos(theta)*sin(phi)*l;
                    ly=sin(theta)*sin(phi)*l;
                    lz=cos(phi)*l;
                else
                    lx=cos(theta)*l;
                    ly=sin(theta)*l;
                    lz=[];
                    [xt,yt,~]=get_coordinates_on_virtual_sample(x+lx,y+ly,[], im_par);
                    ind = sub2ind(sz,round(xt),round(yt)); % the associated index in virtual sample
                    if isempty(intersect(ind,w))  % Is ind on a different pattern ?
                        stick_out=1;
                    end
                    n_trial=n_trial+1;
                end
            end
            if stick_out==1 % Success
                % Keep these values
                sm{lx_idx}=lx;
                sm{ly_idx}=ly;
                sm{lz_idx}=lz;
            else
                SMISMessage1=['Linkage error cannot be set for fluorophore: ', sm_par.fluorophore_name,'  !'];
                SMISMessage2='Increase linkage length or remove ''sticking out'' option';
                SMISMessage=[SMISMessage1 newline SMISMessage2];
                warndlg(SMISMessage,'Warning')
                uiwait
                x2=x;
                y2=x;
                z2=x;
                ok=0;
                return
            end
        else % Just use the fixed linkage error
            lx=sm{lx_idx};
            ly=sm{lx_idx};
            lz=sm{lx_idx};
        end
    else % Case of random linkage error through data collection
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

% Assign the linkage errors
x2=x+lx;
y2=y+ly;
if im_par.simul_3D==1
    z2=z+lz;
else
    z2=[];
end


