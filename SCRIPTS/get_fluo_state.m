function f_out=get_fluo_state(ft_det,state, t1,t2, det_lim)

% PURPOSE:
%   Update detected trace for a certain state, ie in that updated trace, a
%   state is either there during a full frametime, or not
%
% INPUTS:
%	ft_det: the input trace
%	t1: starting time (typically time at start of a frame)
%	t2: end time (typically time at end of a frame)
%   det_lim: [s] minium duration for a state to be detecable
%
% OUTPUTS:
%	f_out: 1 if fluorescent state present during frame, 0 if not.
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2020.

t=ft_det(1,:);
s=ft_det(2,:);

w=find(t>t1 & t<t2);

if isempty(w) % No change during frametime 
    w=find(t<=t1,1,'last');
    if s(w)==state
        f_out=1;
    else
        f_out=0;
    end
    return
else
    w2=find(s(w)==state, 1);
    if isempty(w2) % Fluorescent state does not appear during frametime
        f_out=0;
        return
    else       
        tr_loc=horzcat([t1;s(w(1))],ft_det(:,w),[t2;ft_det(2,w(end))]); % trace between t1 and t2
        
        [state_duration, ~, ~]=get_blink_info(tr_loc(2,:),tr_loc(1,:),state);
        
        if sum(state_duration)>det_lim % Detection
            f_out=1;
        else
            f_out=0;
        end
    end
end



