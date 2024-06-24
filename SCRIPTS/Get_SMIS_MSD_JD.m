function [msd_array,lags_array,jd, mjd,stats] = Get_SMIS_MSD_JD(tracks,frametime,max_lag,pixel_size,mode_3D)

% PURPOSE:
% Produces mean square displacement and jump distances array from Ground truth SMIS simulation track series
%
% INPUTS:
%   tracks: The input tracks
%   frametime: [s] the total frame time (frametime + addtime)
%   max_lag: [s] The maximum lag over which to calculate the MSD's
%   pixel_size: [um] The raster size
%   3D_mode: Set to 1 to perform calculations in 3D
%
% OUTPUTS:
%	msd_array: The array of MSD's, a Nx5 array Containing:
%       1/ the track #;
%       2/ the lags in # of frame;
%       3/ the msds; [um^2]
%       4/ the jump distances (over all lags) [um]
%   lags_array: The array of lag times in [s]
%   jd: The jump distances (array 1x# of lags=1) [um]
%   mjd: The mean jump distances (array 1x# of tracks) [um]
%   stats: Statistics: Structure containing
%       jd_array_mean: Mean jump distance over the different lag times
%       jd_array_std: Standard deviation of jump distance over the different lag times
%       msd_array_mean: Mean MSDs over the different lag times
%       msd_array_std: Standard deviation of MSD's over the different lag times
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2020.
%	D.Bourgeois, January 2022.
%	D.Bourgeois, March 2022. Introduced max_lag to avoid overflow of memory in case of very long tracks !

n_tracks=length(tracks);

if nargin < 5
    mode_3D=0;
end

if n_tracks >0
    
    msd_size=sum([tracks.tracklength].*max_lag); 
    msd_array=zeros(msd_size,4); % 1: track #, 2: Lag time [frames], 3: MSD [um2], 4: Jump distance [um], 5: Mean Jump distance [um]
    msd_id=1; % Initialize msd_id
    
    disp('Extracting MSDs ...')
    MyWaitbar = waitbar(0,'Calculating MSDs ...');
    
    for k=1:n_tracks 
        
        if k/100==fix(k/100)
            waitbar(k/n_tracks,MyWaitbar)
            clc; disp(['Tracks [%]: ',num2str(100*k/n_tracks)]);
        end
        
        u_tr_l=tracks(k).tracklength;
        x=pixel_size*[tracks(k).x]; % convert to microns
        y=pixel_size*[tracks(k).y];
        if mode_3D==1
            z=pixel_size*[tracks(k).z];
        end
        
        t=[tracks(k).frames];
        
        
        for i=1:u_tr_l-1
            xs=circshift(x,[-i,0]);
            ys=circshift(y,[-i,0]);
            if mode_3D==1
                zs=circshift(z,[-i,0]);
            end
            ts=circshift(t,[0, -i]);
            
            lag_t=ts(1:end-i)-t(1:end-i);
            msd_x=((xs(1:end-i)-x(1:end-i)).^2);
            msd_y=((ys(1:end-i)-y(1:end-i)).^2);
            if mode_3D==1
                msd_z=((zs(1:end-i)-z(1:end-i)).^2);
            end
            
            msd_x=msd_x(lag_t<=max_lag);
            msd_y=msd_y(lag_t<=max_lag);
            if mode_3D==1
                msd_z=msd_z(lag_t<=max_lag);
            end
            lag_t=lag_t(lag_t<=max_lag);
            
            n_lags=numel(lag_t);
            
            if mode_3D==0
                msd=msd_x+msd_y;
            else
                msd=msd_x+msd_y+msd_z;
            end
            
            msd_array(msd_id:msd_id+n_lags-1,1)=k; % The track #
            msd_array(msd_id:msd_id+n_lags-1,2)=lag_t; % The lags in # of frame
            msd_array(msd_id:msd_id+n_lags-1,3)=msd; % The msds
            msd_array(msd_id:msd_id+n_lags-1,4)=sqrt(msd); % The jump distances
            
            msd_id=msd_id+n_lags; % Increment msd_id
        end
    end
    
    msd_array=msd_array(msd_array(:,1)>0,:); % Make sure msd_array contains only useful data
    
    close(MyWaitbar)

    %% Calculate the MJD's
    disp('Extracting MJDs ...')
    lags=msd_array(:,2);
    jd=msd_array(lags==1,4);
    jd_tr=msd_array(lags==1,1); % The corresponding track ids
    n_tr_used=numel(unique(jd_tr));
    mjd=zeros(1,n_tr_used);
    for loc_i=1:n_tr_used
        mjd(loc_i)=mean(jd(jd_tr==loc_i));
    end
    mjd=mjd(~isnan(mjd));
    
    
    %% Calculate statistics
    disp('Extracting statistics ...')
    target_lags=frametime*(1:max_lag); % Lags to investigate [s]
    n_lags=numel(target_lags); % # of lags to investigate
    lags_array=frametime*(1:n_lags);
    
    msd_array_mean=zeros(1,n_lags);
    msd_array_std=zeros(1,n_lags);
    jd_array_mean=zeros(1,n_lags);
    jd_array_std=zeros(1,n_lags);
    
    
    for i=1:n_lags
        w=find(abs((msd_array(:,2)*frametime)-target_lags(i)) < 1e-10);
        if ~isempty(w)
            msd_distribution=msd_array(w,3);
            msd_array_mean(i)=mean(msd_distribution);
            msd_array_std(i)=std(msd_distribution)/sqrt(length(msd_distribution)); % standard error of the mean
            jd_distribution=msd_array(w,4);
            jd_array_mean(i)=mean(jd_distribution);
            jd_array_std(i)=std(jd_distribution)/sqrt(length(jd_distribution)); % standard error of the mean
        else
            msd_array_mean(i)=NaN;
            msd_array_std(i)=NaN;
            jd_array_mean(i)=NaN;
            jd_array_std(i)=NaN;
        end
    end
    
    % remove points for which there are no data
    lags_array=lags_array(~isnan(msd_array_mean));
    msd_array_std=msd_array_std(~isnan(msd_array_mean));
    msd_array_mean=msd_array_mean(~isnan(msd_array_mean));
    jd_array_std=jd_array_std(~isnan(jd_array_mean));
    jd_array_mean=jd_array_mean(~isnan(jd_array_mean));
    stats.msd_array_mean=msd_array_mean;
    stats.msd_array_std=msd_array_std;
    stats.jd_array_mean=jd_array_mean;
    stats.jd_array_std=jd_array_std;
    
    % cut the MSD array
    w=msd_array(:,2)<=max_lag;
    msd_array=msd_array(w,:);
else
    msd_array=[];
    lags_array=[];
    stats=[];
end

disp('Done !')

end