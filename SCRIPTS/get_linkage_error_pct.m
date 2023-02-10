function [x2,y2,z2,sm,ok] = get_linkage_error_pct(x,y,z,sm, im_par, sm_pattern_indices, sm_par)
%
% NAME:
%	get_linkage_error
% PURPOSE:
%
% INPUTS:
%   x, y, z: coordinate of the single molecule on detector image [raster units]
%   sm: the single molecule (with coordinates on high-resolution image) in raster units
%	im_par: the imaging parameters
%   sm_pattern_indices: indices of virtual sample subpatterns
%	sm_par: the sm parameters
%
% OUTPUTS:
%	x2,y2,z2: the corrected x y z for linkage error [raster units]
%   sm: The updated single molecule for linkage error
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2022
%	D.Bourgeois, February 2023, introduce sm_pattern_indices, now disconnected from sm_par

% Get the proper indices
c_sp_idx=1;
lx_idx=2;
ly_idx=3;
lz_idx=4;
le_set_idx=5;
ok=1;

if sm_par.linkage_pattern_control==0 % Case of free linkage error
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
else % Case of linkage error forcing the fluorophore to stick to specific pattern

    max_trials=1000; % Maximun number of trials

    if sm_par.linkage_fixed==1 % Case of fixed linkage error through data collection
        if isempty(sm{lx_idx}) % Set linkage error once
            pattern_found=0;
            n_trial=1;
            % Get the target sub pattern
            w_cp_id=find(sm_par.n_sp_id==sm{c_sp_idx});
            w=sm_pattern_indices.w_patterns(sm_par.n_sp_id==sm_par.linkage_pattern_sp_id(w_cp_id)).w; % Indices on the target subpattern
            if im_par.simul_3D==1 % 3D Case
                sz=im_par.binning*[im_par.n, im_par.m, im_par.nz]; % Size of virtual sample image
            else
                sz=im_par.binning*[im_par.n, im_par.m]; % Size of virtual sample image
            end

            l=(sm_par.linkage_length+randn*sm_par.linkage_std)/im_par.raster; % Length of linkage error [raster unit]

            while pattern_found==0 && n_trial<max_trials
                theta=2*pi*rand; % Pick random angle
                if im_par.simul_3D==1 % 3D Case
                    phi=pi*rand; % Pick random angle
                    lx=cos(theta)*sin(phi)*l;
                    ly=sin(theta)*sin(phi)*l;
                    lz=cos(phi)*l;
                    [xt,yt,zt]=get_coordinates_on_virtual_sample(x+lx,y+ly,z+lz, im_par);
                    ind = sub2ind(sz,round(xt),round(yt),round(zt)); % the associated index in virtual sample
                    if ~isempty(find(w==ind,1))  % Is ind on the desired pattern ?
                        pattern_found=1;
                    end
                    n_trial=n_trial+1;
                else % 2D Case
                    lx=cos(theta)*l;
                    ly=sin(theta)*l;
                    lz=[];
                    [xt,yt,~]=get_coordinates_on_virtual_sample(x+lx,y+ly,[], im_par);
                    ind = sub2ind(sz,round(xt),round(yt)); % the associated index in virtual sample
                    if ~isempty(find(w==ind,1))  % Is ind on the desired pattern ?
                        pattern_found=1;
                    end
                    n_trial=n_trial+1;
                end
            end
            if pattern_found==0 % Failure
                SMISMessage1=['Linkage error cannot be set for fluorophore: ', sm_par.fluorophore_name,'  !'];
                SMISMessage2='Increase linkage length or remove ''Control linkage pattern_'' option';
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
        if sm{le_set_idx}==im_par.current_frame % The linkage error has been already calculated for the other channel
            lx=sm{lx_idx};
            ly=sm{lx_idx};
            lz=sm{lx_idx};
        else % Repeat the procedure over
            pattern_found=0;
            n_trial=1;
            % Get the target sub pattern
            w_cp_id=find(sm_par.n_sp_id==sm{c_sp_idx});
            w=sm_pattern_indices.w_patterns(sm_par.n_sp_id==sm_par.linkage_pattern_sp_id(w_cp_id)).w; % Indices on the target subpattern
            if im_par.simul_3D==1 % 3D Case
                sz=im_par.binning*[im_par.n, im_par.m, im_par.nz]; % Size of virtual sample image
            else
                sz=im_par.binning*[im_par.n, im_par.m]; % Size of virtual sample image
            end

            l=(sm_par.linkage_length+randn*sm_par.linkage_std)/im_par.raster; % Length of linkage error [raster unit]

            while pattern_found==0 && n_trial<max_trials
                theta=2*pi*rand; % Pick random angle
                if im_par.simul_3D==1 % 3D Case
                    phi=pi*rand; % Pick random angle
                    lx=cos(theta)*sin(phi)*l;
                    ly=sin(theta)*sin(phi)*l;
                    lz=cos(phi)*l;
                    [xt,yt,zt]=get_coordinates_on_virtual_sample(x+lx,y+ly,z+lz, im_par);
                    ind = sub2ind(sz,round(xt),round(yt),round(zt)); % the associated index in virtual sample
                    if ~isempty(find(w==ind,1))  % Is ind on the desired pattern ?
                        pattern_found=1;
                    end
                    n_trial=n_trial+1;
                else % 2D Case
                    lx=cos(theta)*l;
                    ly=sin(theta)*l;
                    lz=[];
                    [xt,yt,~]=get_coordinates_on_virtual_sample(x+lx,y+ly,[], im_par);
                    ind = sub2ind(sz,round(xt),round(yt)); % the associated index in virtual sample
                    if ~isempty(find(w==ind,1))  % Is ind on the desired pattern ?
                        pattern_found=1;
                    end
                    n_trial=n_trial+1;
                end
            end
            if pattern_found==0 % Failure
                SMISMessage1=['Linkage error cannot be set for fluorophore: ', sm_par.fluorophore_name,'  !'];
                SMISMessage2='Increase linkage length or remove ''Control linkage pattern_'' option';
                SMISMessage=[SMISMessage1 newline SMISMessage2];
                warndlg(SMISMessage,'Warning')
                uiwait
                x2=x;
                y2=x;
                z2=x;
                ok=0;
                return
            end
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



