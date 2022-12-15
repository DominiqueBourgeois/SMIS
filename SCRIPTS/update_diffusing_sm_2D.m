function  sm = update_diffusing_sm_2D(sm, im_par, sm_par, S_DS)

% PURPOSE:
%	Update position and diffusion state of diffusing SM that are in
%	starting diffusion state S_DS
%
% INPUTS:
%   sm: the single molecules (with coordinates on high-resolution image)
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%   S_DS: starting diffusion state
%
% OUTPUTS:
%   sm: the single molecules updated for: diffusion states, current position
%
% MODIFICATION HISTORY:
%	D.Bourgeois, December 2019: version > simulate_palm_vsn15.3
%	D.Bourgeois, November 2020: version > simulate_palm_vsn16.3 (introduce
%	Velocity)
%	D.Bourgeois, January 2022: modify script to account for fast exchange 
%       rates between diffusion states in the diffuse PSF mode. Adapted to SMIS 1.3
%	D.Bourgeois, April 2022: option to record the whole diffusion state
%       history, but not used at this stage.

record_ds_history=0; % Do not record the ds history for SMIS 1.3. May be useful for a future version

% set general parameters
dt_addtime=im_par.dt_diff_addtime;
dt_frametime=im_par.dt_diff_frametime;
dt=dt_addtime+dt_frametime;

% dye_mode=sm_par.mode;
raster=im_par.raster/im_par.binning; % pixel size in high res (true) image [nm]
S=[im_par.binning*im_par.n, im_par.binning*im_par.m]; % size of high-resolution image

% D: diffusion coefficients [nm^2/s]
% C: diffusion confinement
% K: diffusion exchange rate matrix
D=sm_par.D; %[um2.s-1]
V=sm_par.V; %[um.s-1]
C=sm_par.D_confined;
K=sm_par.D_ex_rates;
%Convert D in units of raster: D is given in nm^2/s
D_ras0=1e+06*D/raster^2; % Ex: if D = 0.01 um2/s = 100 nm^2/s = (10 nm x 10 nm)/s, and pixel size = 100 nm, in pixel D = 0.01 pix^2/s
%Convert V in units of raster: V is given in um/s
V_ras0=1e+03*V/raster;

x_h=[sm.x]; % x current position of the molecules
% x_h_d=x_h; % replicate for moved x
y_h=[sm.y]; % y current position of the molecules
% y_h_d=y_h; % replicate for moved y

% Check if exchange rates are compatible with frame time and add time
if im_par.use_diffuse_psf==1
    if ~isempty(sm_par.D_ex_rates)
        max_ex_rate=max(sm_par.D_ex_rates(:,3));
        %Get the maximum subframe time [s] with oversampling of 10
        max_dt=1/(max_ex_rate*im_par.ex_rates_min_oversampling);
    else
        max_dt=inf;
    end
else
    max_dt=[]; % Needed for parallel computing
end

% Needed for parallel computing
current_frame=im_par.current_frame;
use_diffuse_psf=im_par.use_diffuse_psf;
diffuse_psf_radius=im_par.diffuse_psf_radius;

parfor (k=1:numel(sm), im_par.parforArg) % Go for all activated molecules
% for k=1:numel(sm)
    V_ras=V_ras0; % Reassign to avoid broadcasting error in parfor 
    D_ras=D_ras0; % Reassign to avoid broadcasting error in parfor 
    
    %Initialize track if molecule activated in current frame
    if sm(k).activated==current_frame
        sm(k).x_track=[x_h(k),current_frame];
        sm(k).y_track=[y_h(k),current_frame];
    end
    
    % In current position, check the potential evolution of the SM
    
    % Do precise calculations if use_diffuse_psf=1
    if use_diffuse_psf==1
        
        % First let the sm move during addtime
        if dt_addtime > 0
            if max_dt>=dt_addtime % In that case no need to subdivide
                % Look for a potential change in diffusion state
                [sm(k).n_diff_state, sm(k).n_sp]=get_pattern_transition([x_h(k),y_h(k)], ...
                    V_ras(S_DS), D_ras(S_DS), S_DS,K,C,dt_addtime,sm_par,im_par); % Get new diffusion coefficient
                
                % Update position and check if a potential change in diffusion state is realized or not
                [x_h_d_tmp, y_h_d_tmp,sm(k)]=get_new_XY(x_h(k), y_h(k), sm(k), D_ras, dt_addtime, S, sm_par, im_par);

                if record_ds_history==1
                    % Diffusion state history during addtime
                    ds_add_h=[[0,S_DS];[dt_addtime,sm(k).diff_state]];
                end
                
            else % In that case we have to subdivide
                n_substeps=ceil(dt_addtime/max_dt);
                dt_subaddtime=dt_addtime/n_substeps;
                sub_xy=zeros(n_substeps+1,2); % define the intermediate positions sub_x and sub_y
                sub_xy(1,:)=[x_h(k), y_h(k)]; % starting position
                c_diff_state=S_DS; % Current diffusion state
                
                if record_ds_history==1
                    ds_add_h=nan(n_substeps+1,2); % Diffusion state history during addtime
                    ds_add_h(1,:)=[0,S_DS]; % Assign first state to starting state at time 0
                end
                
                for k2=1:n_substeps
                    % Look for a potential change in diffusion state
                    [sm(k).n_diff_state, sm(k).n_sp]=get_pattern_transition([sub_xy(k2,1),sub_xy(k2,2)], ...
                        V_ras(c_diff_state), D_ras(c_diff_state), c_diff_state,K,C,dt_subaddtime,sm_par,im_par); % Get new diffusion coefficient
                    
                    % Update position and check if a potential change in diffusion state is realized or not
                    [sub_xy(k2+1,1),sub_xy(k2+1,2),sm(k)]=get_new_XY(sub_xy(k2,1),sub_xy(k2,2), sm(k), ...
                        D_ras, dt_subaddtime, S, sm_par, im_par);
                    
                    % Reassign current diffusion state
                    c_diff_state=sm(k).diff_state;
                    
                    if record_ds_history==1
                        % Diffusion state history during addtime
                        ds_add_h(k2+1,:)=[k2*dt_subaddtime,c_diff_state];
                    end
                end
                x_h_d_tmp=sub_xy(end,1); % ending position
                y_h_d_tmp=sub_xy(end,2);
            end
            
            N_DS=sm(k).diff_state; % Reassign starting diffusion state
            
        else % Nothing needs to be changed
            x_h_d_tmp=x_h(k);
            y_h_d_tmp=y_h(k);
            N_DS=S_DS; % Do no change the starting diffusion state
            
            if record_ds_history==1
                % Diffusion state history during addtime
                ds_add_h=[];
            end
        end
        
        %Then evaluate how it moves during frametime
        
        % First do a quick evaluation to estimate the diffused distance
        % Look for a potential change in diffusion state
        [sm(k).n_diff_state, sm(k).n_sp]=get_pattern_transition([x_h_d_tmp,y_h_d_tmp], ...
            V_ras(N_DS), D_ras(N_DS), N_DS,K,C,dt_frametime,sm_par,im_par); 
        
        % Update position and check if a potential change in diffusion state is realized or not
        [x_h_d_tmp2, y_h_d_tmp2,sm_tmp2]=get_new_XY(x_h_d_tmp, y_h_d_tmp, ...
            sm(k), D_ras, dt_frametime, S, sm_par, im_par);
        
        diffused_distance=raster*sqrt((x_h_d_tmp2-x_h_d_tmp)^2+(y_h_d_tmp2-y_h_d_tmp)^2); % in [nm]
        
        % Then start sub-calculations if the diffused distance is too big or if exchange rates between diffusion states is too fast
        if diffused_distance>diffuse_psf_radius || max_dt<dt_frametime 
            %# of substeps
            n_substeps=ceil(max([diffused_distance/diffuse_psf_radius,dt_frametime/max_dt]));
            dt_subframetime=dt_frametime/n_substeps;
            sub_xy=zeros(n_substeps+1,2); % define the intermediate positions sub_x and sub_y
            sub_xy(1,:)=[x_h_d_tmp, y_h_d_tmp]; % starting position
            
            c_diff_state=N_DS; % Current diffusion state
            
            if record_ds_history==1
                ds_fr_h=nan(n_substeps+1,2); % Diffusion state history during addtime
                ds_fr_h(1,:)=[dt_addtime,c_diff_state]; % Assign first state to starting state
            end
            
            for k2=1:n_substeps
                % Look for a potential change in diffusion state
                [sm(k).n_diff_state, sm(k).n_sp]=get_pattern_transition([sub_xy(k2,1),sub_xy(k2,2)], ...
                    V_ras(c_diff_state), D_ras(c_diff_state), c_diff_state,K,C,dt_subframetime,sm_par,im_par); % Get new diffusion coefficient
                
                % Update position and check if a potential change in diffusion state is realized or not
                [sub_xy(k2+1,1),sub_xy(k2+1,2),sm(k)]=get_new_XY(sub_xy(k2,1),sub_xy(k2,2), sm(k), ...
                    D_ras, dt_subframetime, S, sm_par, im_par);
                
                % Reassign current diffusion state
                c_diff_state=sm(k).diff_state;
                
                if record_ds_history==1
                    % Diffusion state history during frametime
                    ds_fr_h(k2+1,:)=[dt_addtime+k2*dt_subframetime,c_diff_state];
                end
            end
            x_h_d=sub_xy(end,1); % ending position
            y_h_d=sub_xy(end,2);
        else % if moved distance was small enough, or if there is no issue with fast exchange rates keep it as is
            x_h_d=x_h_d_tmp2; % Keep initial diffused position
            y_h_d=y_h_d_tmp2;
            sm(k)=sm_tmp2; % Keep initial updated sm
            
            if record_ds_history==1
                % Diffusion state history during frametime
                ds_fr_h=[[dt_addtime,N_DS];[dt_addtime+dt_frametime,sm(k).diff_state]];
            end
        end
        
        if record_ds_history==1
            %Get the entire diffusion state history
            if ~isempty(ds_add_h)
                ds_h=vertcat(ds_add_h,ds_fr_h(2:end,:)); % Full diffusion state history during addtime + frametime
            else
                ds_h=ds_fr_h;
            end
        end
        
    else % case of no diffuse PSF
        
        % Look for a change in diffusion state
        [sm(k).n_diff_state, sm(k).n_sp]=get_pattern_transition([x_h(k),y_h(k)], V_ras(S_DS), D_ras(S_DS), S_DS,K,C,dt,sm_par,im_par); % Get new diffusion coefficient

        % Update position and check if a potential change in diffusion state is realized or not
        [x_h_d, y_h_d,sm(k)]=get_new_XY(x_h(k), y_h(k), sm(k), D_ras, dt, S, sm_par, im_par);
        
        if record_ds_history==1
            ds_h=[[0,S_DS];[dt_addtime+dt_frametime,sm(k).diff_state]]; % Full diffusion state history during addtime + frametime
        end
    end
    
    if record_ds_history==1
        % Take diffusion state at half addtime+frametime
        %ds_t2=(ds_h(:,1)-(dt_addtime+dt_frametime)/2).^2;
        %w_midtime=find(ds_t2==min(ds_t2),1);
        
        %if isempty(sm(k).diff_state_trace) % Initialize the diffusion state trace at the beginning
        %    sm(k).diff_state_trace=[current_frame-1;sm(k).init_diff_state];
        %end
        %sm(k).diff_state_trace=horzcat(sm(k).diff_state_trace, [current_frame;ds_h(w_midtime,2)]);
        sm(k).diff_state_trace=vertcat(sm(k).diff_state_trace, ds_h);
    else
        sm(k).diff_state_trace=horzcat(sm(k).diff_state_trace, [current_frame;sm(k).diff_state]);
    end
    
    % Update coordinates on high resolution image
    sm(k).x=x_h_d;
    sm(k).y=y_h_d;
    % Update sub_coordinates on high resolution image if
    % diffuse_psf
    if use_diffuse_psf==1
        if diffused_distance>diffuse_psf_radius
            sm(k).sub_x=sub_xy(:,1);
            sm(k).sub_y=sub_xy(:,2);
        else
            sm(k).sub_x=[];
            sm(k).sub_y=[];
        end
    end
    % Update diffusion track for next frame
    sm(k).x_track=vertcat(sm(k).x_track, [x_h_d,current_frame+1]);
    sm(k).y_track=vertcat(sm(k).y_track, [y_h_d,current_frame+1]);
end

