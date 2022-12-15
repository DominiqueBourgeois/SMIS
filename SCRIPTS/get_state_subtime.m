function [t_on, fr] =get_state_subtime(trace, state, n_subtime)

%
% PURPOSE:
%   Get the times a trace is in a particular state, for each n_subtime
%   segment (Case of a diffuse PSF)
%
% INPUTS:
%   trace: the single molecule trace: array(2,N) with array(1,:)=the time
%   in [s] and array(2,:)=the corresponding states
%   state: the state to be looked at
%   n_subtime: the number of equally spaced segments in which the trace will be divided
%
% OUTPUTS:
%   fr: the fraction of total time spent in specified state for each segment
%   t_on: total time spent in specified state
%
% MODIFICATION HISTORY:
%	D.Bourgeois, October 2019.

% Look if that state is found in the trace
if ~isempty(trace)
    w=find(trace(2,:)==state);
    % If yes, proceed
    if ~isempty(w)
        sub_t=zeros(n_subtime,1);
        t=trace(1,w);
        t_Start=t(1:2:end); % All times where that state starts
        t_End=t(2:2:end); % All times where that state stops
        T=trace(1,end); % Lengths of the trace

        % Go over all segments
        for i=1:n_subtime
            t1=(i-1)*T/n_subtime; % Time start for the segment
            t2=i*T/n_subtime; % Time end for the segment

            % Shrink t_Start to relevant times for the segment
            ts=t_Start;
            ts(t_Start<=t1 | t_Start>t2)=[];

            % Shrink t_End to relevant times for the segment
            te=t_End;
            te(t_End<=t1 | t_End>t2)=[];

            %Get state at start of segment
            seg_start=[]; % Set seg_start

            %Last start of state
            ls=t1-t_Start; ls=t_Start(ls>=0);
            if ~isempty(ls)
                ls=ls(end);
            else
                seg_start=0; % if there was no start before beginning of segment, then the state is off at the beginning
            end

            if isempty(seg_start) % If there was a start before beginning of segment, look whether it ended or not before beginning of segment
                %Last end of state
                le=t1-t_End; le=t_End(le>=0);
                if ~isempty(le) % If the state ended before beginning of segment, We need to analyze
                    le=le(end);
                    if le>ls; seg_start=0; else seg_start=1; end % Look at which event is the latest before beginning of segment
                else seg_start=1; % If the state did not end before beginning of segment, Then the state is on at the beginning
                end
            end

            if seg_start==1 % On time Started before beginning of segment
                if ~isempty(ts) ; ts=horzcat(t1,ts); else ts=t1; end
            end
            if numel(ts)>numel(te) % A on_time finished after end of segment
                if ~isempty(te); te=horzcat(te,t2); else te=t2; end
            end

            % Now this is the on time for the relevant state within the segment
            sub_t(i)=sum(te-ts);
        end
        t_on=sum(sub_t);
        fr=sub_t/t_on;
    else
        t_on=0;
        fr=0;
    end
else
    t_on=0;
    fr=0;
end

