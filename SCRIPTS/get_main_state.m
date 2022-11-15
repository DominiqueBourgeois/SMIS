function ft_det=get_main_state(ft_det,state, t1,t2, det_lim,addtime)

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
%	ft_det: the updated trace
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2020.

t=ft_det(1,:);
s=ft_det(2,:);

w=find(t>t1 & t<t2);

if isempty(w)
    return
else
    w2=find(s(w)==state, 1);
    if isempty(w2)
        ft_det(:,w)=[];
        return
    else

        tr_loc=ft_det(:,w); % trace between t1 and t2
        
        [state_duration, ~, ~]=get_blink_info(tr_loc(2,:),tr_loc(1,:),state);
        
        if sum(state_duration)>det_lim % Detection
            if w(1)>1 && w(end) < numel(s) % There are events before and after the current frame
                if s(w(1)-1)==state && s(w(end)+1)==state % No change
                    ft_det(:,w)=[];
                    return
                end
                if s(w(1)-1)==state && s(w(end)+1)~=state % Change at the end
                    ft_det=horzcat(ft_det(:,1:w(1)-1),[t2;state],[t2;s(w(end)+1)],ft_det(:,w(end)+1:end));
                    return
                end
                if s(w(1)-1)~=state && s(w(end)+1)==state % Change at the beginning
                    ft_det=horzcat(ft_det(:,1:w(1)-1),[t1;s(w(1)-1)],[t1;state],ft_det(:,w(end)+1:end));
                    return
                end
                if s(w(1)-1)~=state && s(w(end)+1)~=state % Change at the beginning and the end
                    ft_det=horzcat(ft_det(:,1:w(1)-1),[t1;s(w(1)-1)],[t1;state],[t2;state],[t2;s(w(end)+1)],ft_det(:,w(end)+1:end));
                    return
                end
            elseif w(1)==1 && w(end) < numel(s) % There are events only after the current frame
                if s(w(end)+1)==state % No change
                    ft_det=horzcat([t1;0], [t1;state], ft_det(:,w(end)+1:end));
                    return
                end
                if s(w(end)+1)~=state % Change at the end
                    ft_det=horzcat([t1;0],[t1;state],[t2;state],[t2;s(w(end)+1)],ft_det(:,w(end)+1:end));
                    return
                end
            elseif w(1)>1 && w(end) == numel(s) % There are events only before the current frame
                if s(w(1)-1)==state  % No change
                    ft_det=horzcat(ft_det(:,1:w(1)-1),[t2;state],[t2;0]);
                    return
                end
                if s(w(1)-1)~=state % Change at the beginning
                    ft_det=horzcat(ft_det(:,1:w(1)-1),[t1;s(w(1)-1)],[t1;state],[t2;state],[t2;0]);
                    return
                end
            elseif  w(1)==1 && w(end) == numel(s) % There are no events before or after the current frame
                ft_det=horzcat([t1;0],[t1;state],[t2;state],[t2;0]);
                return
            end
        else % No detection
            if w(1)>1 && w(end) < numel(s) % There are events before and after the current frame
                if s(w(1)-1)~=state && s(w(end)+1)~=state % No change
                        ft_det(:,w)=[];
                        return
                end
                if s(w(1)-1)~=state && s(w(end)+1)==state % Change at the end, Set the state id during frame to -1
                    ft_det=horzcat(ft_det(:,1:w(1)-1),[t2;0], [t2;state],ft_det(:,w(end)+1:end));
                    return
                end
                if s(w(1)-1)==state && s(w(end)+1)~=state % Change at the beginning Set the state id during frame to -1
                    ft_det=horzcat(ft_det(:,1:w(1)-1),[t1;state], [t1;0],ft_det(:,w(end)+1:end));
                    return
                end
                if s(w(1)-1)==state && s(w(end)+1)==state % Change at the beginning and the end
                    ft_det=horzcat(ft_det(:,1:w(1)-1),[t1;state],[t1;0],[t2;0],[t2;state],ft_det(:,w(end)+1:end));
                    return
                end
            elseif w(1)==1 && w(end) < numel(s) % There are events only after the current frame
                if s(w(end)+1)~=state % No change
                    ft_det(:,w)=[];
                    return
                end
                if s(w(end)+1)==state % Change at the end
                    ft_det=horzcat([t2;0],[t2;state],ft_det(:,w(end)+1:end));
                    return
                end
            elseif w(1)>1 && w(end) == numel(s) % There are events only before the current frame
                if s(w(1)-1)~=state  % No change
                    ft_det(:,w)=[];
                    return
                end
                if s(w(1)-1)==state % Change at the beginning
                    ft_det=horzcat(ft_det(:,1:w(1)-1), [t1;state],[t1;0]);
                    return
                end
            elseif  w(1)==1 && w(end) == numel(s) % There are no events before or after the current frame
                ft_det=[];
                return
            end
        end
    end
    ft_det=tr_loc;
end



