function split_tracks = SplitTracks(tracks,t,simul_3D)

% PURPOSE:
% Split tracks into shorter tracks for means of getting as many as possible short but continuous short tracks
% of constant length
%
% INPUTS:
%   tracks: The input tracks
%   t: The desired tracklength
%   simul_3D: Set to one if a 3D track
%
% OUTPUTS:
%   split_tracks: The output split tracks
%
% MODIFICATION HISTORY:
%	D.Bourgeois, January 2022.


n_tracks=length(tracks);

if n_tracks >0
    
    %Initialize the short tracks
    split_tracks=struct;
    
    
    id=1;
    % Go over all tracks
    for k=1:n_tracks
        
        % Split the track into continuous segments
        [on_times, ~, ~, start_frames]=get_on_off_times_from_frames_on(tracks(k).frames,0);
        
        %Number of segments
        n_seg=numel(on_times);
        for j=1:n_seg
            
            % Get the local continuous track
            w=find(tracks(k).frames==start_frames(j),1);
            loc_track.x=tracks(k).x(w:w+on_times(j)-1);
            loc_track.y=tracks(k).y(w:w+on_times(j)-1);
            if simul_3D==1
                loc_track.z=tracks(k).z(w:w+on_times(j)-1);
            end
            loc_track.frames=tracks(k).frames(w:w+on_times(j)-1);
            loc_track.tracklength=on_times(j);
            
            % Now split this local continues tracks into shorter tracks
            tl=loc_track.tracklength;
            
            nst=fix((tl-1)/t); % # of split tracks
            if nst>0
                for i=1:nst
                    si=1+(t-1)*(i-1); % Start index
                    ei=1+(t-1)*(i-1)+t-1; % End index
                    split_tracks(id).x=loc_track.x(si:ei);
                    split_tracks(id).y=loc_track.y(si:ei);
                    if simul_3D==1
                        split_tracks(id).z=loc_track.z(si:ei);
                    end
                    split_tracks(id).frames=loc_track.frames(si:ei);
                    split_tracks(id).tracklength=t;
                    split_tracks(id).ref_track_id=k;
                    id=id+1;
                end
            end
        end
    end
else
    split_tracks=[];
end

end