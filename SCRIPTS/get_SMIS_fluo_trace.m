function frames_on=get_SMIS_fluo_trace(ft,fluo_states, frametime,addtime, n_images, min_detectable_on)

% PURPOSE:
%   Get detected fluorescence trace for a certain state
%
% INPUTS:
%	ft: initial fluorescence trace
%	fluo_states: the fluorescent state(s) id(s)
%   frametime: frametime in [s]
%   n_images: the number of images in the dataset
%   addtime: addtime in [s]
%   min_detectable_on: minimum fraction of time the state should be present during frametime to be considered detectable
%
% OUTPUTS:
%	frames_on: list of frames during which dye is in specified state
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2020.
%	D.Bourgeois, January 2022.



if size(ft,2)==1 && any(nnz(~(ft(2,1)-fluo_states))) % Deal with the case where a molecule does not change state until the end
    ft=[ft,ft]; % Add a point at the end of the dataset
    ft(1,end)=(frametime+addtime-1e-9)*n_images; % Remove 1e-9 to avoid f_end > n_images
end

t_start=ft(1,1); % Start of trace [s]
t_end=ft(1,end);

f_start=1+fix(t_start/(frametime+addtime)); % First frame
f_end=1+fix(t_end/(frametime+addtime));

frames_on=zeros(2,f_end-f_start+1);
frames_on(1,:)=f_start:f_end;

acqu_time=n_images*(frametime+addtime);

% Go through all frames during which state appears
for k=f_start:f_end
    t1=(k-1)*(frametime+addtime)+addtime; % Only look during frametime, not during addtime
    t2=k*(frametime+addtime);
    frames_on(2,k-f_start+1)=get_SMIS_fluo_state(ft,fluo_states, t1,t2,min_detectable_on*frametime,acqu_time);
end


