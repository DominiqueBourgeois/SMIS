function [tracks, subtracks, off_times]=extract_SMIS_GT_tracks(all_sm, par)

% NAME:
%	extract_ground_truth_tracks
%
% PURPOSE:
%   Get the ground truth track from a sm
% INPUTS:
%   all_sm: The single molecules, possibly merged from several fluorophore types if this is a multicolor experiment
%   fluorescent_states: The ids of the fluorescent state producing the signal
%   par: Input parameters
%       blinking_gap: Maximum # of dark frames in tracks
%       plot_tracks: plot [plot_tracks]% of the tracks
%       color_track_mode: 1: single random color per track ; 0: color gradient per track
%       safety_margin: % [pixels] safety margin for molecules going out of the FOV (this can create artificial blinks: use -999 to not supress parts of the tracks out of the FOV).
%       min_phot: % Minimum number of photons for detection
%       select_bleached: % Set to 1 to only select bleached molecules
%       select_activated: % Set to 1 to only select activated molecules
%       im_par: % The imaging parameters
%       Extractsubtracks: % Set to 1 to extract subtracks as well
%
% OUTPUTS:
%   tracks: All the identified GT tracks, there might be more tracks than molecules if molecules go out of the field of view and return into it.
%   subtracks: All the identified subtracks due to blinking: off time in a track > blinking_gap will split a track into two sub tracks.
%   off_times: All the identified of times [frames]
%
% MODIFICATION HISTORY:
%	D.Bourgeois, January 2019.
%	D.Bourgeois, March 2020. Compatibility with PALM simulation software vsn15.2
%	D.Bourgeois, October 2020. Compatibility with PALM simulation software vsn16.2
%	D.Bourgeois, January 2022. Compatibility with SMIS 1.3, include segments for different diffusion coefficients
%	D.Bourgeois, March 2022. Remove segments extraction (done instead in extract_SMIS_segments.m) and option to extract subtracks
%	D.Bourgeois, October 2023. %Treat the eventual case where x_track not recorded in a frame_on event (ie photons detected before photoconversion)

im_par=par.im_par;
% plot_tracks=par.plot_tracks;
% color_track_mode=par.color_track_mode;
% plot_tracks=0;
% color_track_mode=0;
safety_margin=par.safety_margin;
ch=par.channel;

if isempty(all_sm)
    disp('No molecule was selected !')
    tracks=[];
    subtracks=[];
    off_times=[];
    return
end

% Get the # of mol for each molecule type
total_n_mol=size(all_sm,2);

off_times=zeros(1,10*total_n_mol);

tracks(1:10*total_n_mol)=struct(... % There might be more tracks than total number of molecules if molecules transiently go out of the field of view.
    'x',[], ...
    'y',[], ...
    'z',[], ...
    'frames',[], ...
    'tracklength', 0, ...
    'tracklength_with_blink', 0, ...
    'n_subtracks', 0, ...
    'subtracks_frame_start', [], ...
    'subtracks_frame_end', [], ...
    'subtracks_length', [], ...
    'off_times', [], ...
    'segments', [], ... % Parts of the track corresponding to different diffusion states
    'molecule_id', [], ... % The id of the molecule forming the track
    'fluorophore_id', [] ... % For multicolor experiments
    );

if par.Extractsubtracks==1
    subtracks(1:50*total_n_mol)=struct(...
        'x',[], ...
        'y',[], ...
        'z',[], ...
        'frames',[], ...
        'tracklength_with_blink', [], ... % # of frames from beginning to end of track
        'tracklength', [], ...          % # of measured points
        'off_times', [], ...
        'segments', [], ... % Parts of the subtrack corresponding to different diffusion states
        'molecule_id', [], ... % The id of the molecule forming the track
        'fluorophore_id', [] ... % For multicolor experiments
        );
else
    subtracks=[];
end

k_o=1; % Counter for off_times
k2=1; % index for subtracks
k3=1; % index for subtracks
k_fast_bleach=0; % index for # of molecules that bleach too fast to produce track

% if plot_tracks>0; figure; hold on; end

n_within_FOV=0;
n_escape_FOV=0;
%evaluate tracks
disp('Extracting tracks ...')
MyWaitbar = waitbar(0,'Calculating tracks ...');

for k=1:total_n_mol % Go over total number of molecules
    
    if k/100==fix(k/100)
        waitbar(k/total_n_mol,MyWaitbar)
        clc; disp(['Tracks [%]: ',num2str(100*k/total_n_mol)]);
    end
    
    sm=all_sm(k);
    
    molecule_id=k;
    
    if ch==1
        fluo_frames=sm.frames_on_ch1(sm.n_phot_det_ch1>par.min_phot);
    else
        fluo_frames=sm.frames_on_ch2(sm.n_phot_det_ch2>par.min_phot);
    end
    
    %look if mol has been activated
    if ~isempty(fluo_frames)
        if numel(fluo_frames)>1 % only proceed if fluorescence observed in at least 2 frames
            frames_on=fluo_frames;
            
            % Get x and y for visible part of track
            if isempty(sm.x_track) % If track is empty, activation & bleaching occured in the same frame
                sm.x_track=[sm.x,sm.bleached];
                sm.y_track=[sm.y,sm.bleached];
                if im_par.simul_3D==1
                    sm.z_track=[sm.z,sm.bleached];
                end
            end
            
            if max(sm.x_track(:,2))> fluo_frames(end) % If track exceed max fluo_frame, the sm probably did not bleach, remove last point in track
                sm.x_track(end,:)=[];
                sm.y_track(end,:)=[];
                if im_par.simul_3D==1
                    sm.z_track(end,:)=[];
                end
            end
            
            x=sm.x_track(:,1);
            y=sm.y_track(:,1);
            if im_par.simul_3D==1
                z=sm.z_track(:,1);
            end
            track_f=sm.x_track(:,2); % the tracks frames
            [~, ix]=intersect(track_f,frames_on);

            %Treat the eventual case where x_track not recorded in a
            %frame_on event (ie photons detected before photoconversion)
            if numel(ix)~=numel(frames_on)
                % Frames unique to frames_on
                [missed_frames,missed_frames_ix]=setdiff(frames_on,track_f);
                disp(['Molecule: ',num2str(k),' fluoresced at frames: ',num2str(missed_frames),' before track started to be recorded !']);
                % Let's ignore those frames
                frames_on(missed_frames_ix)=[];
            end
            
            %Get coordinates on detector (x and y are in high-res image from
            %simulations
            x_frames=(x(ix)-0.5)/im_par.binning+0.5; % x coordinate in EMCCD image
            y_frames=(y(ix)-0.5)/im_par.binning+0.5; % y coordinate in EMCCD image
            if im_par.simul_3D==1
                z_frames=(z(ix)-0.5)/im_par.binning+0.5; % z coordinate in EMCCD image;
            end
            
            % Eventually plot track
            %             if plot_tracks>0
            %                 diary off
            %                 if k*plot_tracks/100==fix(k*plot_tracks/100)
            %                     clc; disp(['Tracks [%] : ',num2str(100*k/total_n_mol)]);
            %                     if color_track_mode==1 % single color per track
            %                         ColorVector = rand(1,3);
            %                         plot(x_frames_red_on, y_frames_red_on, 'Color', ColorVector, 'LineWidth',2);
            %                     elseif color_track_mode==0 %  color gradient per track
            %                         ColorVector = ColorCoding(full_tracklength);
            %                         for n = 1 : full_tracklength-1
            %                             plot(x_frames_red_on(n:n+1), y_frames_red_on(n:n+1), 'Color', ColorVector(n, :), 'LineWidth',2);
            %                         end
            %                     end
            %                     title('Plot of full ground truth tracks')
            %                 end
            %                 diary off
            %             end
            
            % Look if mol did get out of the FOV, in which case molecule is not excited by laser anymore => nonsense
            % then select parts of the track that are within FOV
            
            %case when mol stays within FOV
            if min(x_frames)>=(0.5+safety_margin) && max(x_frames)<=(im_par.n+0.5-safety_margin) && min(y_frames)>=(0.5+safety_margin) && max(y_frames)<=(im_par.m+0.5-safety_margin)
                
                n_tracks_within_FOV=1;
                loc_tracks=struct('x',[],'y',[],'z',[],'frames',[],'tracklength', 0);
                loc_tracks.x=x_frames;
                loc_tracks.y=y_frames;
                if im_par.simul_3D==1
                    loc_tracks.z=z_frames;
                end
                loc_tracks.frames=frames_on;
                loc_tracks.tracklength=length(frames_on);
                loc_tracks.tracklength_with_blink=max(frames_on)-min(frames_on)+1;
                loc_tracks.molecule_id=molecule_id;
                n_within_FOV=n_within_FOV+1;
                FOV_check=1; % Flag to say this is OK
            else % the mol transiently or irreversibely goes out of FOV
                w_in_FOV=find(x_frames>=(0.5+safety_margin) & x_frames<=(im_par.n+0.5-safety_margin) & y_frames>=(0.5+safety_margin) & y_frames<=(im_par.m+0.5-safety_margin));
                if ~isempty(w_in_FOV) % in case molecule never goes into the FOV (this should never happen)
                    x_in_FOV=x_frames(w_in_FOV);
                    y_in_FOV=y_frames(w_in_FOV);
                    if im_par.simul_3D==1
                        z_in_FOV=z_frames(w_in_FOV);
                    end
                    frames_in_FOV=frames_on(w_in_FOV);
                    
                    shifted_w_in_FOV=circshift(w_in_FOV,[0,-1])-w_in_FOV;
                    shifted_w_in_FOV=shifted_w_in_FOV(1:end-1);
                    w_new_track_within_FOV=find(shifted_w_in_FOV>1);
                    n_escape_FOV=n_escape_FOV+1;
                    if ~isempty(w_new_track_within_FOV) % transient escape from FOV
                        
                        n_tracks_within_FOV=length(w_new_track_within_FOV)+1; % # of parts within the FOV
                        loc_tracks(1:n_tracks_within_FOV)=struct(...
                            'x',[],...
                            'y',[],...
                            'z',[],...
                            'frames',[],...
                            'tracklength', 0,...
                            'fluorophore_id', [],...
                            'molecule_id', []...
                            );
                        
                        %subtrack #1
                        loc_tracks(1).x=x_in_FOV(1:w_new_track_within_FOV(1));
                        loc_tracks(1).y=y_in_FOV(1:w_new_track_within_FOV(1));
                        if im_par.simul_3D==1
                            loc_tracks(1).z=z_in_FOV(1:w_new_track_within_FOV(1));
                        end
                        loc_tracks(1).frames=frames_in_FOV(1:w_new_track_within_FOV(1));
                        loc_tracks(1).tracklength=length(loc_tracks(1).frames);
                        loc_tracks(1).tracklength_with_blink=max([loc_tracks(1).frames])-min([loc_tracks(1).frames])+1;
                        loc_tracks(1).molecule_id=molecule_id;
                        if n_tracks_within_FOV>2
                            for loc_i=1:length(w_new_track_within_FOV)-1
                                loc_tracks(loc_i+1).x=x_in_FOV(w_new_track_within_FOV(loc_i)+1:w_new_track_within_FOV(loc_i+1));
                                loc_tracks(loc_i+1).y=y_in_FOV(w_new_track_within_FOV(loc_i)+1:w_new_track_within_FOV(loc_i+1));
                                if im_par.simul_3D==1
                                    loc_tracks(loc_i+1).z=z_in_FOV(w_new_track_within_FOV(loc_i)+1:w_new_track_within_FOV(loc_i+1));
                                end
                                loc_tracks(loc_i+1).frames=frames_in_FOV(w_new_track_within_FOV(loc_i)+1:w_new_track_within_FOV(loc_i+1));
                                loc_tracks(loc_i+1).tracklength=length(loc_tracks(loc_i+1).frames);
                                loc_tracks(loc_i+1).tracklength_with_blink=max([loc_tracks(loc_i+1).frames])-min([loc_tracks(loc_i+1).frames])+1;
                                loc_tracks(loc_i+1).molecule_id=molecule_id;
                            end
                        end
                        %last subtrack
                        loc_tracks(end).x=x_in_FOV(w_new_track_within_FOV(end)+1:end);
                        loc_tracks(end).y=y_in_FOV(w_new_track_within_FOV(end)+1:end);
                        if im_par.simul_3D==1
                            loc_tracks(end).z=z_in_FOV(w_new_track_within_FOV(end)+1:end);
                        end
                        loc_tracks(end).frames=frames_in_FOV(w_new_track_within_FOV(end)+1:end);
                        loc_tracks(end).tracklength=length(loc_tracks(end).frames);
                        loc_tracks(end).tracklength_with_blink=max([loc_tracks(end).frames])-min([loc_tracks(end).frames])+1;
                        loc_tracks(end).molecule_id=molecule_id;
                    else % irreversible escape
                        n_tracks_within_FOV=1; % # of parts within the FOV
                        loc_tracks=struct('x',[],'y',[],'z',[],'frames',[],'tracklength', 0);
                        loc_tracks(1).x=x_in_FOV;
                        loc_tracks(1).y=y_in_FOV;
                        if im_par.simul_3D==1
                            loc_tracks(1).z=z_in_FOV;
                        end
                        loc_tracks(1).frames=frames_in_FOV;
                        loc_tracks(1).tracklength=length(loc_tracks(1).frames);
                        loc_tracks(1).tracklength_with_blink=max([loc_tracks(1).frames])-min([loc_tracks(1).frames])+1;
                        loc_tracks(1).molecule_id=molecule_id;
                    end
                    FOV_check=1;
                else
                    FOV_check=0;
                    n_escape_FOV=n_escape_FOV+1;
                end
            end
            
            if  FOV_check==1
                for loc_i=1:n_tracks_within_FOV
                    tracks(k2).x=[loc_tracks(loc_i).x];
                    tracks(k2).y=[loc_tracks(loc_i).y];
                    if im_par.simul_3D==1
                        tracks(k2).z=[loc_tracks(loc_i).z];
                    end
                    tracks(k2).frames=[loc_tracks(loc_i).frames];
                    tracks(k2).tracklength=loc_tracks(loc_i).tracklength;
                    tracks(k2).tracklength_with_blink=max([tracks(k2).frames])-min([tracks(k2).frames])+1;
                    tracks(k2).molecule_id=molecule_id;
                    
                    frames_on=[loc_tracks(loc_i).frames];
                    x_frames=[loc_tracks(loc_i).x];
                    y_frames=[loc_tracks(loc_i).y];
                    if im_par.simul_3D==1
                        z_frames=[loc_tracks(loc_i).z];
                    end
                    
                    frame_red_on_shifted=circshift(frames_on,[0,-1]);
                    delta_frames=frame_red_on_shifted(1:end-1)-frames_on(1:end-1)-1;
                    
                    w_blink=find(delta_frames>=1); % these are the  off times > max_blink
                    
                    if ~isempty(w_blink)
                        off_time_frame_start=frames_on(w_blink);
                        off_time_frame_end=frames_on(w_blink+1);
                        blinks=off_time_frame_end-off_time_frame_start-1;
                        off_times(k_o:k_o+length(w_blink)-1)=blinks;
                        tracks(k2).off_times=blinks;
                        k_o=k_o+length(w_blink); % update Counter for off_times
                    end
                    
                    %Now get the subtracks
                    if par.Extractsubtracks==1
                        
                        
                        w_long_delta=find(delta_frames>par.blinking_gap); % these are the  off times > max_blink
                        
                        if ~isempty(w_long_delta)
                            tracks(k2).n_subtracks=length(w_long_delta)+1;
                            subtracks_frame_start=[frames_on(1,1),frames_on(w_long_delta+1)];
                            subtracks_frame_end=[frames_on(w_long_delta),frames_on(1,end)];
                            subtracks_lengths=subtracks_frame_end-subtracks_frame_start+1;
                            tracks(k2).subtracks_frame_start=subtracks_frame_start;
                            tracks(k2).subtracks_frame_end=subtracks_frame_end;
                            tracks(k2).subtracks_length=subtracks_lengths; % this includes the possible dark frames
                            %tracks(k2).fluorophore_id=fluorophore_id;
                            tracks(k2).molecule_id=molecule_id;
                            for i=1:tracks(k2).n_subtracks
                                if i==1
                                    ind=1:w_long_delta(1);
                                elseif i<tracks(k2).n_subtracks
                                    ind=(w_long_delta(i-1)+1):w_long_delta(i);
                                else
                                    ind=(w_long_delta(i-1)+1):tracks(k2).tracklength;
                                end
                                subtracks(k3).x=x_frames(ind);
                                subtracks(k3).y=y_frames(ind);
                                if im_par.simul_3D==1
                                    subtracks(k3).z=z_frames(ind);
                                end
                                subtracks(k3).frames=frames_on(ind);
                                subtracks(k3).tracklength_with_blink=subtracks_lengths(i);
                                subtracks(k3).tracklength=length(subtracks(k3).frames);
                                subtracks(k3).molecule_id=molecule_id;
                                sub_frames_on=frames_on(ind);
                                sub_frame_red_on_shifted=circshift(sub_frames_on,[0,-1]);
                                sub_delta_frames=sub_frame_red_on_shifted(1:end-1)-sub_frames_on(1:end-1)-1;
                                
                                w_sub_blink=find(sub_delta_frames>=1); % these are the  off times in subtracks
                                if ~isempty(w_sub_blink)
                                    sub_off_time_frame_start=sub_frames_on(w_sub_blink);
                                    sub_off_time_frame_end=sub_frames_on(w_sub_blink+1);
                                    sub_blinks=sub_off_time_frame_end-sub_off_time_frame_start-1;
                                    subtracks(k3).off_times=sub_blinks;
                                end
                                k3=k3+1;
                            end
                        else
                            subtracks(k3).x=tracks(k2).x;
                            subtracks(k3).y=tracks(k2).y;
                            if im_par.simul_3D==1
                                subtracks(k3).z=tracks(k2).z;
                            end
                            subtracks(k3).frames=tracks(k2).frames;
                            subtracks(k3).tracklength_with_blink=max([tracks(k2).frames])-min([tracks(k2).frames])+1;
                            subtracks(k3).tracklength=tracks(k2).tracklength;
                            subtracks(k3).molecule_id=molecule_id;
                            subtracks(k3).off_times=tracks(k2).off_times;
                            k3=k3+1;
                        end
                    end
                    k2=k2+1;
                end
            end
        else
            k_fast_bleach=k_fast_bleach+1;
        end
    end
end

close(MyWaitbar)
disp(['Number of tracks simulated: ', num2str(total_n_mol)])
disp(['Number of tracks that stayed within FOV: ', num2str(n_within_FOV)])
disp(['Number of tracks that escaped FOV: ', num2str(n_escape_FOV)])
disp(['Number of molecules that bleached too fast to produce track: ', num2str(k_fast_bleach)])

off_times(k_o:end)=[];
tracks([tracks.tracklength]==0)=[]; % remove non existing tracks

if par.Extractsubtracks==1
    subtracks(k3:end)=[]; % supress non existing subtracks
    subtracks([subtracks.tracklength]==0)=[]; % remove non existing subtracks
end

disp(['Total number of blinks: ', num2str(k_o-1)])
disp(['Number of individual tracks output (some simulated tracks go in and out of FOV several times: ', num2str(length(tracks))])
disp(['Number of individual subtracks output (due to intermittencies > max_blink): ', num2str(length(subtracks))])
disp(['Number of tracks: ', num2str(numel(tracks))])
disp(['Number of molecules: ', num2str(total_n_mol)]);


% if plot_tracks>0; hold off; end
end