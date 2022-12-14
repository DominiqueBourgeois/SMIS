function sm=get_fluorescence_pct(sm, sm_par, im_par)

%
% PURPOSE:
%   Get the fluorescence (on_times, emitted counts, emitted spectra, of a sm
%   for all fluorescent states during current frame
%
% INPUTS:
%   sm: the single molecule in cell format
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   sm: the single molecules updated for fluorescence
%
% MODIFICATION HISTORY:
%	D.Bourgeois, August 2019.
%	D.Bourgeois, September 2022. Adapted for pct toolbox

%Here are the indices in sm cell array
sub_x_idx = 1;
c_sp_idx = 2;
fluo_trace_idx = 3;
Ns_idx = 4;
sampling_rate_idx = 5;
n_abs_idx = 6;
n_em_idx = 7;
tot_n_em_idx = 8;
bleached_idx = 9;
t_on_idx = 10;
fr_t_on_idx = 11;

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
            sm{t_on_idx}(i)=get_state_time_pct(sm{fluo_trace_idx},sm{sampling_rate_idx},sm_par.fluorescent_states(i));
        else % Case for a diffuse PSF in case of diffusion
            n_subtime=numel(sm{sub_x_idx});
            [sm{t_on_idx}(i),sm{fr_t_on_idx}(i,1:n_subtime)] =get_state_subtime(sm{fluo_trace_idx},sm_par.fluorescent_states(i),n_subtime);
        end
        % Get # of absorbed photons during that time
        n_abs=sm{t_on_idx}(i)*sm{Ns_idx}(end,sm_par.fluorescent_states(i))*sm{sampling_rate_idx}(end);
        sm{n_abs_idx}(sm_par.fluorescent_states(i))=n_abs;
        
        % Get # of emitted photons during that time; use poissonian statistics
        if n_abs>0
            if sm_par.fluorogenic==1
                %n_em=poissrnd(sm_par.quantum_yield(i)*n_abs*sm_par.fluorogenicity(sm_par.n_sp_id==sm{c_sp_idx}));
                n_em=sm_par.quantum_yield(i)*n_abs*sm_par.fluorogenicity(sm_par.n_sp_id==sm{c_sp_idx}); % Poissonian process will be applied later
            else
                %n_em=poissrnd(sm_par.quantum_yield(i)*n_abs);
                n_em=sm_par.quantum_yield(i)*n_abs; % Poissonian process will be applied later
            end
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






