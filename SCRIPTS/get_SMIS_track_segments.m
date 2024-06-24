function segments=get_SMIS_track_segments(track,sm,n_diff_states,par)

% PURPOSE:
%   Get segments in a sm track corresponding to the parts of tracks with defined diffusion state
%
% INPUTS:
%	track: the track
%	sm: the sm that produced the track
%   n_diff_states: the # of diffusion states
%   par: parameters (used for 3D)
%
% OUTPUTS:
%	segments: The different segments [structure]
%
% MODIFICATION HISTORY:
%	D.Bourgeois, January 2022.

% Define segments
segments(1:n_diff_states)=struct('n',[],'on_times',[],'frames',[],'x',[],'y',[],'z',[],'diff_state_id',[],'diff_state_origin',[],'tr_id',[]);

% All diffusion states and frames present in track (even during off times)
all_frames=sm.diff_state_trace(1,:);
all_ds=sm.diff_state_trace(2,:);

% Types of diffusion states that were present in the track
ds_id=unique(all_ds);

% Go through each diffusion state
for k=1:numel(ds_id)
    frames_ds=all_frames(all_ds==ds_id(k)); % Frames where the ds is there
    
    [~, ~, ids]=intersect(frames_ds,track.frames); % Ids in track of those frames that were detected with current ds   
    frames_on_ds=track.frames(ids); % The detected frames for that particular diffusion state
    x_on=track.x(ids);
    y_on=track.y(ids);
    if par.im_par.simul_3D==1
        z_on=track.z(ids);
    end
    
    % See how many segments we have
    [seg_on_times, ~,~,start_frames]=get_on_off_times_from_frames_on(frames_ds,0); % Use a gap of 0
    
    n_seg=numel(seg_on_times); % Number of segments with that diffusion state over the whole track, observed or not
    % Define structure for segments
    % Be very careful to clear the structures first
    clear('frames','x','y','z');
    frames(1:n_seg)=struct('values',[]); % Detected frames for each segment
    x(1:n_seg)=struct('values',[]); % x-values for each segment
    y(1:n_seg)=struct('values',[]); % y-values for each segment
    if par.im_par.simul_3D==1
        z(1:n_seg)=struct('values',[]); % z-values for each segment
    end
    obs_on_times=nan(1,n_seg); % The observed on times for each segment
    orig=nan(1,n_seg); % Origin diffusion state for each transition
       
    
    for i=1:n_seg
        %Detected frames corresponding to the on time
        [frames(i).values, sub_ids]=intersect(frames_on_ds, start_frames(i):(start_frames(i)+seg_on_times(i)-1));
   
        if ~isempty(sub_ids) % Only fill the values if some frames were on during this diffusion state lifetime
            x(i).values=x_on(sub_ids);
            y(i).values=y_on(sub_ids);
            if par.im_par.simul_3D==1
                z(i).values=z_on(sub_ids);
            end
            
            % Find out the observed on time
            obs_on_times(i)=frames(i).values(end)-frames(i).values(1)+1;
        end
        % Find out the origin state
        if start_frames(i)>1
            orig_state=all_ds(all_frames==start_frames(i)-1);
            if ~isempty(orig_state) % orig_state might be empty if start_frames(i)=sm.activated
                orig(i)=orig_state;
            end
        end
    end
     
    % Cleanup
    w_empty=find(isnan(obs_on_times));
    frames(w_empty)=[]; %#ok<*AGROW>
    x(w_empty)=[];
    y(w_empty)=[];
    if par.im_par.simul_3D==1
        z(w_empty)=[];
    end
    obs_on_times(w_empty)=[];
    obs_orig=orig;
    obs_orig(w_empty)=[];
    
    n_ds_on=numel(frames); % Actual number of observed segments
        
    segments(ds_id(k)).diff_state_id=ds_id(k); % Id of the diffusion state producing those segments
    segments(ds_id(k)).n=n_ds_on; % Number of observed segments with that diffusion state
    segments(ds_id(k)).n_tot=n_seg; % Number of segments with that diffusion state, observed or not
    segments(ds_id(k)).on_times_tot=seg_on_times; % On times for the segments in [frames], observed or not 
    segments(ds_id(k)).on_times=obs_on_times; % On times for the observed segments in [frames]
    segments(ds_id(k)).frames=frames; % Frames in the observed segment
    segments(ds_id(k)).x=x; % x-values in the observed segment
    segments(ds_id(k)).y=y; % y-values in the observed segment
    if par.im_par.simul_3D==1
        segments(ds_id(k)).z=z; % z-values in the observed segment
    end
    segments(ds_id(k)).diff_state_origin_tot=orig; % Diffusion state that just precedes,observed or not 
    segments(ds_id(k)).diff_state_origin=obs_orig; % Diffusion state that just precedes,observed 
    segments(ds_id(k)).tr_id=par.k; % Keep in memory the ID of the track from which the segment originates
end

%Treat particular case of transition during first frame
% if sm.init_diff_state~=all_ds(1)
%     segments(sm.init_diff_state).diff_state_id=sm.init_diff_state; % Id of the diffusion state producing those segments
%     segments(sm.init_diff_state).tr_id=par.k; % Keep in memory the ID of the track from which the segment originates
%     
%     if ~isempty(segments(sm.init_diff_state).n) % There are segments with this ds
%         segments(sm.init_diff_state).n=segments(sm.init_diff_state).n+1; % Increase by 1
%         segments(sm.init_diff_state).n_tot=segments(sm.init_diff_state).n_tot+1; % Increase by 1
%         segments(sm.init_diff_state).on_times_tot=[segments(sm.init_diff_state).on_times_tot, 0]; % Add an ontime of zero
%         segments(sm.init_diff_state).on_times=[segments(sm.init_diff_state).on_times, 0]; % Add an ontime of zero
%         
%         frames=segments(sm.init_diff_state).frames;
%         frames(numel(frames)+1).values=[]; % Add empty frame
%         segments(sm.init_diff_state).frames=frames; % Frames in the observed segment
%         
%         x=segments(sm.init_diff_state).x;
%         x(numel(x)+1).values=[]; % Add empty x
%         segments(sm.init_diff_state).x=x; % x-values in the observed segment
%         
%         y=segments(sm.init_diff_state).y;
%         y(numel(y)+1).values=[]; % Add empty y
%         segments(sm.init_diff_state).y=y; % y-values in the observed segment
%                         
%         if par.im_par.simul_3D==1
%             z=segments(sm.init_diff_state).z;
%             z(numel(z)+1).values=[]; % Add empty z
%             segments(sm.init_diff_state).z=z; % z-values in the observed segment
%         end
%         
%         segments(sm.init_diff_state).diff_state_origin_tot=[segments(sm.init_diff_state).diff_state_origin_tot, nan]; % Diffusion state that just precedes,observed or not
%         segments(sm.init_diff_state).diff_state_origin=[segments(sm.init_diff_state).diff_state_origin, nan]; % Diffusion state that just precedes,observed
%         
%     else
%         segments(sm.init_diff_state).n=1; % Increase by 1
%         segments(sm.init_diff_state).n_tot=1; % Increase by 1
%         segments(sm.init_diff_state).on_times_tot=0; % Add an ontime of zero
%         segments(sm.init_diff_state).on_times=0; % Add an ontime of zero
%         
%         segments(sm.init_diff_state).frames.values=[]; % Frames in the observed segment
%         segments(sm.init_diff_state).x.values=[]; % x-values in the observed segment
%         
%         segments(sm.init_diff_state).y.values=[]; % y-values in the observed segment                
%         if par.im_par.simul_3D==1
%             segments(sm.init_diff_state).z.values=[]; % z-values in the observed segment
%         end
%         
%         segments(sm.init_diff_state).diff_state_origin_tot=nan; % Diffusion state that just precedes,observed or not
%         segments(sm.init_diff_state).diff_state_origin=nan; % Diffusion state that just precedes,observed
%     end
%         
% end



