function ft_det=get_detected_fluo_trace2(ft,state, frametime,addtime, min_detectable_on)

% PURPOSE:
%   Get detected fluorescence trace for a certain state
%
% INPUTS:
%	ft: initial fluorescence trace
%	state: the state id
%   frametime: frametime in [s]
%   addtime: addtime in [s]
%   min_detectable_on: minimum fraction of time the state should be present during frametime to be considered detectable
%
% OUTPUTS:
%	ft_det: the detected trace for the state in question
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2020.


ft_det=ft;

t_start=ft(1,1); % Start of trace [s]
t_end=ft(1,end);

tot_frametime=frametime+addtime;

f_start=1+fix(t_start/tot_frametime); % First frame
f_end=1+fix(t_end/tot_frametime);

n_frames=f_end-f_start+1;
% frames_on=zeros(2,n_frames);
% frames_on(1,:)=f_start:f_end;

T=[]; % Initialize T trace

% Go through all frames during which state appears: this gets a list of
% frames during which molecule is considered on
for k=f_start:f_end
    t1=(k-1)*tot_frametime+addtime;
    t2=k*tot_frametime;
    f_out=get_fluo_state(ft,state, t1,t2,min_detectable_on*frametime);
    T=horzcat(T,[(k-1)*tot_frametime;f_out],[t2;f_out]);
end

%transform T into an output trace with correct format.
% for k=1:n_frames
TS=circshift(T,[0,1]);
T2=T-TS;
w_on=find(T2(2,:)==1);
w_off=find(T2(2,:)==-1);
w_on(w_on==1)=[];
w_off(w_off==1)=[];

ft_det=[]; % Initialize with starting state
if ~isempty(w_on)
    for k=1:numel(w_on)
        ft_det=horzcat(ft_det,T(:,w_on(k)-1:w_on(k)));
    end
end

if ~isempty(w_off)
    for k=1:numel(w_off)
        ft_det=horzcat(ft_det,T(:,w_off(k)-1:w_off(k)));
    end
end

if ~isempty(ft_det)
    [~,si]=sort(ft_det(1,:));
    ft_det=ft_det(:,si);
end

ft_det=horzcat(T(:,1), ft_det, T(:,end));


