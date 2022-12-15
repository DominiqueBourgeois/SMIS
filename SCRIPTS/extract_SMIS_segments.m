function [segments, tracks]=extract_SMIS_segments(all_sm,tracks,par)

% PURPOSE:
%   Extract and merge all segments from a set of tracks
%
% INPUTS:
%   all_sm: The single molecules, possibly merged from several fluorophore types if this is a multicolor experiment
%	tracks: the set of tracks
%   par: Input parameters
%
% OUTPUTS:
%	segments: The merged segments [structure]
%	tracks: The updated tracks [structure]
%
% MODIFICATION HISTORY:
%	D.Bourgeois, January 2022.
%	D.Bourgeois, March 2022.

if par.n_diff_states>1 % Only do it if necessary, ie there are selected molecule with multiple diffusion states

    MyWaitbar = waitbar(0,['Extracting segments for ',par.title,' ...']);

    for k=1:numel(tracks)
        if k/100==fix(k/100)
            waitbar(k/numel(tracks),MyWaitbar)
            clc; disp(['Fraction of ',par.title,' [%]: ',num2str(100*k/numel(tracks))]);
        end
        molecule_id=tracks(k).molecule_id;
        sm=all_sm(molecule_id);
        par.k=k; % Keep track of the track ID
        tracks(k).segments=get_SMIS_track_segments(tracks(k),sm,par.n_diff_states,par);
    end
    close(MyWaitbar) 
    
    
%     %Then do it for the subtracks
%     disp('Extracting segments for subtracks ...')
%     MyWaitbar = waitbar(0,'Extracting segments for subtracks  ...');
%     for k=1:numel(subtracks)
%         if k/100==fix(k/100)
%             waitbar(k/numel(subtracks),MyWaitbar)
%             clc; disp(['Tracks [%]: ',num2str(100*k/numel(subtracks))]);
%         end
%         
%         molecule_id=subtracks(k).molecule_id;
%         sm=all_sm(molecule_id);
%         par.k=k; % Keep track of the track ID
%         subtracks(k).segments=get_SMIS_track_segments(subtracks(k),sm,n_diff_states,par);
%     end
%     close(MyWaitbar)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%



% Get the total number of segments
n_seg=0;
for i=1:par.n_diff_states
    MyData = arrayfun(@(s)s.segments(i).n,tracks,'UniformOutput',false);
    n=[MyData{:}];
    n_seg=n_seg+sum(n);
end

%Define the segments, they will contain the same field as tracks, so that they can be processed as tracks later on
segments(1:n_seg)=struct();

%Go through all tracks
seg_id=1; % Initialize the segment ID number
for i=1:par.n_diff_states % Go through all diffusion states
    for j=1:numel(tracks)
        s=tracks(j).segments(i);
        if ~isempty(s.n)
            for k=1:s.n
                if s.on_times(k)==0
                    disp('h')
                end
                segments(seg_id).x=s.x(k).values;
                segments(seg_id).y=s.y(k).values;
                if ~isempty(s.z)
                    segments(seg_id).z=s.z(k).values;
                else
                    segments(seg_id).z=[];
                end
                segments(seg_id).frames=s.frames(k).values;
                segments(seg_id).tracklength=length(s.frames(k).values);
                if ~isempty(segments(seg_id).frames)
                    segments(seg_id).tracklength_with_blink=max([segments(seg_id).frames])-min([segments(seg_id).frames])+1;
                    [segments(seg_id).on_times,segments(seg_id).off_times,~]=get_on_off_times_from_frames_on(segments(seg_id).frames,0);
                else
                    segments(seg_id).tracklength_with_blink=0;
                    segments(seg_id).on_times=0;
                    segments(seg_id).off_times=[];
                end
                                    
                segments(seg_id).diff_state_id=s.diff_state_id;     
                segments(seg_id).diff_state_origin=s.diff_state_origin(k);
                segments(seg_id).tr_id=s.tr_id;
                segments(seg_id).molecule_id=tracks(j).molecule_id;
                segments(seg_id).fluorophore_id=tracks(j).fluorophore_id;
                
                seg_id=seg_id+1; % Increment the segment ID number
                
            end
        end
    end
end


