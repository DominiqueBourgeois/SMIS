function [state_duration, state_origin, state_occurence, active_time]=get_state_info(s,t,state_id,acq_time, merge)

% PURPOSE:
%   Get info on state duration and origin state for a certain photophysical
%   state (or set of states) with id state_id.
%
% INPUTS:
%	s: state trace from sm
%	t: time trace from sm
%   state_id: state(s) id for the state
%   acq_time: total acquisition time of data set [s]
%   merge: set to 1 to merge different states and consider them as a unique super-state"
%
% OUTPUTS:
%	state_duration: duration of the state [s]
%   state_origin: state id at the origin of the state
%   state_occurence: # of occurence of the state
%   active_time: total time the state stays, including possible
%   interruptions
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2020.
%	D.Bourgeois, September 2020. Simplified routine
%   D.Bourgeois, December 2021. Designed from the get_blink_info.m routine
%   D.Bourgeois, March 2023. Add active_time

[~, ui]=unique(t);
s=s(ui);
t=t(ui);
if nargin==5
 t_save=t;
end


% The following is to handle on or off times that are a combination of different states
if numel(state_id)==1
    w=find(s==state_id); % Lifetime of state under consideration is typically t(w)-t(w-1)
elseif numel(state_id)==2
    w=find(s==state_id(1) | s==state_id(2)); 
elseif numel(state_id)==3
    w=find(s==state_id(1) | s==state_id(2) | s==state_id(3)); 
elseif numel(state_id)==4
    w=find(s==state_id(1) | s==state_id(2) | s==state_id(3) | s==state_id(4)); 
elseif numel(state_id)==5
    w=find(s==state_id(1) | s==state_id(2) | s==state_id(3) | s==state_id(4) | s==state_id(5)); 
elseif numel(state_id)==6
    w=find(s==state_id(1) | s==state_id(2) | s==state_id(3) | s==state_id(4) | s==state_id(5) | s==state_id(6)); 
elseif numel(state_id)==7
    w=find(s==state_id(1) | s==state_id(2) | s==state_id(3) | s==state_id(4) | s==state_id(5) | s==state_id(6) | s==state_id(7));
elseif numel(state_id)==8
    w=find(s==state_id(1) | s==state_id(2) | s==state_id(3) | s==state_id(4) | s==state_id(5) | s==state_id(6) | s==state_id(7) | s==state_id(8)); 
else
    error('Can not handle more than 8 states ...');
end
  
%Treat the case of several merge states to obtain global on times (or offtimes)
state_no_change=nan;
if nargin==5
    if merge==1
        if ~isempty(w)
            s(w)=state_id(1); % Assign all the positions in the state trace having one of the state specified by rs to the first state rs(1)
            cs=s-circshift(s,[0,-1]); % Remove all the positions in between starting and ending global state except the first one
            cs(1)=nan;
            t(cs==0)=[];
            s(cs==0)=[];
            w=find(s==state_id(1)); % Now look again for the global state
            if numel(s)==1
                state_no_change=1;
            end
        end
    end
end

if ~isempty(w)
    if w(1)>1 % Normal case, the state appears at some point but is not present at start
        state_origin=s(w-1);
        state_duration=t(w)-t(w-1);
        active_time=t(w(end))-t(w(1));
    elseif numel(w)>1 % The first occurence of state does not correspond to a blink
        state_origin=s(w(2:end)-1);
        state_duration=t(w(2:end))-t(w(2:end)-1);
        active_time=t(w(end))-t(w(1));
    elseif numel(w)==1 % The state starts but does not end (no bleaching) and origin is undefined
        state_origin=[];
        state_duration=acq_time-t(1);
        active_time=acq_time-t(1);
    elseif state_no_change==1
        state_duration=t_save(end)-t_save(1); % The whole trace
        state_origin=nan;
        active_time=t_save(end)-t_save(1);
    else
        state_origin=[];
        state_duration=[];
        active_time=[];
    end
    
    state_origin(state_origin==state_id(1))=nan; % Treat the case of the starting state in the trace: the origin of the state transition cannot be itself
    state_occurence=numel(state_duration);
else
        state_duration=[];
        state_origin=[];
        state_occurence=0;
        active_time=[];
end
