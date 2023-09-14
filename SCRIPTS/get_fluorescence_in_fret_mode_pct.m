function sm=get_fluorescence_in_fret_mode_pct(sm, sm_par, im_par)

%
% PURPOSE:
%   Get the fluorescence (on_times, emitted counts, emitted spectra, of a sm
%   for all fluorescent states during current frame, taking into account
%   FRET
%
% INPUTS:
%   sm: the single molecule 
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   sm: the single molecules updated for fluorescence
%
% MODIFICATION HISTORY:
%	D.Bourgeois, August 2019.
%	D.Bourgeois, July 2020.
%	D.Bourgeois, September 2022. Adapted for pct toolbox
%	D.Bourgeois, September 2023. Corrected bug line 68: sm{sampling_rate_idx}(1) 
% instead of sm{sampling_rate_idx}(end) which was reffering to addtime sampling rate 

%Here are the indices in sm cell array
sub_x_idx = 1;
fluo_trace_idx = 2;
Ns_idx = 3;
sampling_rate_idx = 4;
n_abs_idx = 5;
n_em_idx = 6;
tot_n_em_idx = 7;
bleached_idx = 8;
t_on_idx = 9;
fr_t_on_idx = 10;
given_fret_photons_idx=11;
received_fret_photons_idx=12;
matched_idx=13;

plot_fluo_traces=0; % Set to 1 to see fluorescent state appears

% Case for a diffuse PSF in case of diffusion: Reset these variables
if im_par.add_diffusion==1 && im_par.use_diffuse_psf==1 && ~isempty(sm{sub_x_idx})
    sm{fr_t_on_idx}=[];
end


% if no photon was absorbed (the molecule was bleached before current frametime): nothing to do
if sm{bleached_idx}~=0 && sm{bleached_idx}<im_par.current_frame
    sm{t_on_idx}=[];
    sm{n_abs_idx}=[];
    sm{n_em_idx}(:)=0;
    return
end

% First get the on-times for the various fluorescent states during current
% frame

for i=1:sm_par.n_fluorescent_states
    if sm_par.fluorescent_states(i)>=sm_par.initial_fluo_state % Only treat accessible fluorescent states
        % Get the time fluorescence is on
        if im_par.add_diffusion==0 || im_par.use_diffuse_psf==0 || isempty(sm{sub_x_idx}) % Normal case
            sm{t_on_idx}(i)=get_state_time_pct(sm{fluo_trace_idx},sm{sampling_rate_idx}(1),sm_par.fluorescent_states(i));
        else % Case for a diffuse PSF in case of diffusion
            n_subtime=numel(sm{sub_x_idx});
            [sm{t_on_idx}(i),sm{fr_t_on_idx}(i,1:n_subtime)] =get_state_subtime(sm{fluo_trace_idx},sm_par.fluorescent_states(i),n_subtime);
        end
        % Get # of absorbed photons during that time
        n_abs=sm{t_on_idx}(i)*sm{Ns_idx}(end,sm_par.fluorescent_states(i))*sm{sampling_rate_idx}(1);
             
        if sm{matched_idx}>0 % The molecule is in interaction with a FRET partner           
            sm{n_abs_idx}(sm_par.fluorescent_states(i))=n_abs-sm{given_fret_photons_idx}(i)+sm{received_fret_photons_idx}(i);
        else
            sm{n_abs_idx}(sm_par.fluorescent_states(i))=n_abs;
        end
               
        % Get # of emitted photons during that time; 
        if sm{n_abs_idx}(sm_par.fluorescent_states(i))>0
            %n_em=poissrnd(sm_par.quantum_yield(i)*sm{n_abs_idx}(sm_par.fluorescent_states(i))); 
            n_em=sm_par.quantum_yield(i)*sm{n_abs_idx}(sm_par.fluorescent_states(i)); % Poissonian process will be applied later
        else 
            n_em=0; 
        end
        sm{n_em_idx}(i)=n_em;
        %Update total # of emitted photons
        sm{tot_n_em_idx}(i)=sm{tot_n_em_idx}(i)+n_em;
    end
end

if plot_fluo_traces==1
    figure(1)
    p=plot(sm{fluo_trace_idx}(1,:),sm{fluo_trace_idx}(2,:));
    set(p,'LineWidth',1.5);
    set(p,'Color','blue');
    ylim([0 max(sm{fluo_trace_idx}(2,:))+1]);
    xlabel('Time [s]','fontsize',8,'fontweight','b')
    ylabel('State [AU]','fontsize',8,'fontweight','b')
    title(gca,'Single Molecule Trace','FontWeight','bold');
    disp(['# of emitted photons in each fluorescent state:' , num2str(sm{n_em_idx})]);
    input('Ok ? ','s');
end






