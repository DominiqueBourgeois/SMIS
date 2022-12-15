function f_out=get_SMIS_fluo_state(ft_det,fluo_states, t1,t2, det_lim, acqu_time)

% PURPOSE:
%   Update detected trace for a certain fluorescent state(s), ie in that updated trace, a
%   state is either there during a full frametime, or not
%
% INPUTS:
%	ft_det: the input trace
%   fluo_states: the fluorescent state(s) to look at
%	t1: starting time (typically time at start of a frame)
%	t2: end time (typically time at end of a frame)
%   det_lim: [s] minium duration for a state to be detecable
%   acqu_time: the acquisition time of the dataset [s]
% OUTPUTS:
%	f_out: 1 if fluorescent state(s) are present during frame, 0 if not.
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2020.

t=ft_det(1,:);
s=ft_det(2,:);

w=find(t>t1 & t<t2);

if isempty(w) % No change during frametime
    w=find(t<=t1,1,'last');
    if any(nnz(~(s(w)-fluo_states)))
        f_out=1;
    else
        f_out=0;
    end
    return
else
    
    
    if numel(fluo_states)==1
        w2=find(s(w)==fluo_states, 1);
    elseif numel(fluo_states)==2
        w2=find(s(w)==fluo_states(1) | s(w)==fluo_states(2), 1);
    elseif numel(fluo_states)==3
        w2=find(s(w)==fluo_states(1) | s(w)==fluo_states(2) | s(w)==fluo_states(3), 1);
    elseif numel(fluo_states)==4
        w2=find(s(w)==fluo_states(1) | s(w)==fluo_states(2) | s(w)==fluo_states(3) | s(w)==fluo_states(4), 1);
    elseif numel(fluo_states)==5
        w2=find(s(w)==fluo_states(1) | s(w)==fluo_states(2) | s(w)==fluo_states(3) | s(w)==fluo_states(4) | s(w)==fluo_states(5), 1);
    else
        error('Can not handle more than 8 states ...');
    end
    
    if isempty(w2) % Fluorescent state does not appear during frametime
        f_out=0;
        return
    else
        tr_loc=horzcat([t1;s(w(1))],ft_det(:,w),[t2;ft_det(2,w(end))]); % trace between t1 and t2
        
        [state_duration, ~, ~]=get_state_info(tr_loc(2,:),tr_loc(1,:),fluo_states,acqu_time);
        
        if sum(state_duration)>det_lim % Detection
            f_out=1;
        else
            f_out=0;
        end
    end
end



