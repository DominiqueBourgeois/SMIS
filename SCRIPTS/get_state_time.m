function t=get_state_time(sm, state)

%
% PURPOSE:
%   Get the time a trace is in a particular state
%
% INPUTS:
%   trace: the single molecule trace: array(2,N) with array(1,:)=the time
%   in [s] and array(2,:)=the corresponding states
%   state: the state to be looked at
%
% OUTPUTS:
%   t: the time spent in specified state
%
% MODIFICATION HISTORY:
%	D.Bourgeois, august 2019.

if ~isempty(sm.fluo_trace)
    w=find(sm.fluo_trace(2,:)==state);
    if ~isempty(w)
        if numel(w)<2
            % This can happen if there is a state change at the end of a
            % subtrace followed by another state change at the beginning of the next
            % subtrace. The trace at the end of first subtrace looks like:
            % state1-state2 and then for the next subtrace it looks like
            % state2-state3. So it goes like state1-state2-state2-state3, but
            % merge_trace removes the last point of the subtrace in the
            % first subtrace, so it goes state1-state2-state3 and state2 appears
            % only once. Then it lasts for one sample time = 1/sm.sampling_rate
            % disp(['Molecule ', num2str(sm.id),': Trace for photophysical state: ', num2str(sm.state), ' appears for only 1 time point']);
            % disp('There was 2 state changes at the end of one frame (last sample) and at the beginning of the next frame (first sample)')
            t=1/sm.sampling_rate(end);
            return
        end
        sw_w=w-circshift(w,[1,1]);
        w2=find(sw_w==1);
        t=sum(sm.fluo_trace(1,w(w2))-sm.fluo_trace(1,w(w2)-1));
    else
        t=0;
    end
else
    t=0;
end
