function tot_state_trace=update_tot_trace(sm, im_par)
%
% INPUTS:
%   sm: the single molecule (with coordinates on high-resolution image) in raster units
%	im_par: the imaging parameters

% MODIFICATION HISTORY:
%	D.Bourgeois, May 2017.
%	D.Bourgeois, July 2019.

state_trace=sm.state_trace;

if im_par.during_frametime==1
    %We add the total time since beginning of experiment + eventual addtime before the frame time
    T=1e-3*(im_par.frametime+im_par.addtime)*(im_par.current_frame-1)+1e-3*im_par.addtime;
    state_trace(1,:)=T+sm.state_trace(1,:);
else
    %We add the total time since beginning of experiment
    T=1e-3*(im_par.frametime+im_par.addtime)*(im_par.current_frame-1);
    state_trace(1,:)=T+sm.state_trace(1,:);
end

if size(state_trace,2)==2 % No state change
    tot_state_trace=sm.tot_state_trace;
elseif size(state_trace,2)>2
    tot_state_trace=horzcat(sm.tot_state_trace,state_trace(:,2:end-1));
end

end




