function  [on_times, off_times, N, start_frames]=get_on_off_times_from_frames_on(frames_on,gap)

% PURPOSE:
%   Extract on times and off times from the list of frames that were on during single molecule acquisition
% INPUTS:
% frames_on: The list of frame numbers
% gap:   The allowed gap [frames] between two frames to reconstitute on times
%
% OUTPUTS:
%   on_times: Array of the on times  [number of frames]
%   off_times: Array of the off times  [number of frames]
%   N: The number of blinks
%   start_frames: The first frame at the beginning of every detected on time
%
% MODIFICATION HISTORY:
%	D.Bourgeois, December 2021.
%	D.Bourgeois, January 2021, Add start_frames variable

if isempty(frames_on)
    on_times=[];
    off_times=[];
    N=[];
    start_frames=[];
    return
end
d=diff(frames_on);
d(d<=(gap+1))=1;
w=find(d>1);
if ~isempty(w)
    on_times=zeros(1,numel(w)+1);
    start_frames=zeros(1,numel(w)+1);
    off_times=frames_on(w+1)-frames_on(w)-1;
    on_times(1)=frames_on(w(1))-frames_on(1)+1;
    start_frames(1)=frames_on(1);
    if numel(w)>1
        for i=2:numel(w)
            on_times(i)=frames_on(w(i))-frames_on(w(i-1)+1)+1;
            start_frames(i)=frames_on(w(i-1)+1);
        end
    end
    on_times(end)=frames_on(end)-frames_on(w(end)+1)+1;
    start_frames(end)=frames_on(w(end)+1);
    N=numel(off_times);
else
    off_times=[];
    on_times=numel(frames_on);
    start_frames=frames_on(1);
    N=0;
end


