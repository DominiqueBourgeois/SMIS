function sm=get_fluorescence(sm, sm_par, im_par)

%
% PURPOSE:
%   Get the fluorescence (on_times, emitted counts, emitted spectra, of a sm
%   for all fluorescent states during current frame
%
% INPUTS:
%   sm: the single molecule (with coordinates on high-resolution image) in raster units
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   sm: the single molecules updated for fluorescence
%
% MODIFICATION HISTORY:
%	D.Bourgeois, August 2019.

plot_fluo_traces=0; % Set to 1 to see fluorescent state appears

% Case for a diffuse PSF in case of diffusion: Reset these variables
if im_par.add_diffusion==1 && im_par.use_diffuse_psf==1 && ~isempty(sm.sub_x)
    sm.fr_t_on=[];
end


% if no photon was absorbed (the molecule was bleached before current frametime): nothing to do
if sm.bleached~=0 && sm.bleached<im_par.current_frame
    sm.t_on=[];
    sm.n_abs=[];
    sm.n_em(:)=0;
    return
end

% First get the on-times for the various fluorescent states during current
% frame

for i=1:sm_par.n_fluorescent_states
    if sm_par.fluorescent_states(i)>=sm_par.initial_fluo_state % Only treat accessible fluorescent states      
        % Get the time fluorescence is on
        if im_par.add_diffusion==0 || im_par.use_diffuse_psf==0 || isempty(sm.sub_x) % Normal case           
            sm.t_on(i)=get_state_time(sm,sm_par.fluorescent_states(i));
        else % Case for a diffuse PSF in case of diffusion
            n_subtime=numel(sm.sub_x);
            [sm.t_on(i),sm.fr_t_on(i,1:n_subtime)] =get_state_subtime(sm.fluo_trace,sm_par.fluorescent_states(i),n_subtime);
        end
        % Get # of absorbed photons during that time
        n_abs=sm.t_on(i)*sm.Ns(end,sm_par.fluorescent_states(i))*sm.sampling_rate(end);
        sm.n_abs(sm_par.fluorescent_states(i))=n_abs;
        
        % Get # of emitted photons during that time; use poissonian statistics
        if n_abs>0
            if sm_par.fluorogenic==1
                n_em=poissrnd(sm_par.quantum_yield(i)*n_abs*sm_par.fluorogenicity(sm_par.n_sp_id==sm.c_sp));
            else
                n_em=poissrnd(sm_par.quantum_yield(i)*n_abs); 
            end
        else
            n_em=0;
        end
        sm.n_em(i)=n_em;

        %Update total # of emitted photons
        sm.tot_n_em(i)=sm.tot_n_em(i)+n_em;
    end
end

if plot_fluo_traces==1
    figure(1)
    p=plot(sm.fluo_trace(1,:),sm.fluo_trace(2,:));
    set(p,'LineWidth',1.5);
    set(p,'Color','blue');
    ylim([0 max(sm.fluo_trace(2,:))+1]);
    xlabel('Time [s]','fontsize',8,'fontweight','b')
    ylabel('State [AU]','fontsize',8,'fontweight','b')
    title(gca,'Single Molecule Trace','FontWeight','bold');
    disp(['# of emitted photons in each fluorescent state:' , num2str(sm.n_em)]);
    input('Ok ? ','s');
end






