function  [on_times, off_times]=get_SMIS_on_times_from_frames_on(frames_on,gap)

% PURPOSE:
%   Extract on times and off times from the list of frames that were on during single molecule acquisition
% INPUTS:
% frames_on: The list of frame numbers
% gap:   The allowed gap between two frames to reconstitute on times
%
% OUTPUTS:
%   on_times: Array of the on times  [number of frames]
%   off_times: Array of the off times  [number of frames]
%
% MODIFICATION HISTORY:
%	D.Bourgeois, December 2021.

d=diff(frames_on);
d(d<=(gap+1))=1;
w=find(d>1);
if ~isempty(w)
    on_times=zeros(1,numel(w)+1);
    off_times=frames_on(w+1)-frames_on(w)-1;
    on_times(1)=frames_on(w(1))-frames_on(1)+1;
    if numel(w)>1
        for i=2:numel(w)
            on_times(i)=frames_on(w(i))-frames_on(w(i-1)+1)+1;
        end
    end
    on_times(end)=frames_on(end)-frames_on(w(end)+1)+1;
else
    off_times=[];
    on_times=numel(frames_on);
end


