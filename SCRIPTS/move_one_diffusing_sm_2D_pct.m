function  sm = move_one_diffusing_sm_2D_pct(sm, sm_par, im_par, S_DS)

% PURPOSE:
%	Update 2D position and diffusion state of a single sm that are is
%	starting diffusion state S_DS
%
% INPUTS:
%   sm: the single molecule
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%   S_DS: starting diffusion state
%
% OUTPUTS:
%   sm: the single molecule updated for: diffusion state, current position
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2022, optimized for parallel computing


%     w_idx = find(matches(sm_par(i).sm_fn,{'x','y','x_track','y_track',...
%         'sub_x','sub_y','v_x','v_y','c_sp', 'n_sp',...
%         'id','bleached','activated','diff_state','n_diff_state',...
%         'diff_state_trace','matched'})==1);

%Extract the useful indices
x_idx=1;
y_idx=2;
x_track_idx=3;
y_track_idx=4;
sub_x_idx=5;
sub_y_idx=6;
% v_x_idx=7;
% v_y_idx=8;
% c_sp_idx=9;
n_sp_idx=10;
% id_idx=11;
% bleached_idx=12;
activated_idx=13;
diff_state_idx=14;
n_diff_state_idx=15;
diff_state_trace_idx=16;
% matched_idx=17;

% set general parameters
dt_addtime=im_par.dt_diff_addtime;
dt_frametime=im_par.dt_diff_frametime;
dt=dt_addtime+dt_frametime;

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
D_ras=1e+06*D/raster^2; % Ex: if D = 0.01 um2/s = 100 nm^2/s = (10 nm x 10 nm)/s, and pixel size = 100 nm, in pixel D = 0.01 pix^2/s
%Convert V in units of raster: V is given in um/s
V_ras=1e+03*V/raster;

x_h=sm{x_idx}; % x current position of the molecule
y_h=sm{y_idx}; % y current position of the molecule

max_dt=sm_par.max_dt; % Maximum subframe time [s]

%Other stuff
current_frame=im_par.current_frame;
use_diffuse_psf=im_par.use_diffuse_psf;
diffuse_psf_radius=im_par.diffuse_psf_radius;
record_ds_history=sm_par.record_ds_history;

%Initialize track if molecule activated in current frame
if sm{activated_idx}==current_frame
    sm{x_track_idx}=[x_h,current_frame];
    sm{y_track_idx}=[y_h,current_frame];
end

% In current position, check the potential evolution of the SM

%% Do precise calculations if use_diffuse_psf=1
if use_diffuse_psf==1
    % First let the sm move during addtime
    if dt_addtime > 0
        if max_dt>=dt_addtime % In that case no need to subdivide
            % Look for a potential change in diffusion state
            [sm{n_diff_state_idx}, sm{n_sp_idx}]=get_pattern_transition([x_h,y_h], ...
                V_ras(S_DS), D_ras(S_DS), S_DS,K,C,dt_addtime,sm_par,im_par); % Get new diffusion coefficient

            % Update position and check if a potential change in diffusion state is realized or not
            [x_h_d_tmp, y_h_d_tmp,sm]=get_new_XY_pct(x_h, y_h, sm, D_ras, dt_addtime, S, sm_par, im_par);

            if record_ds_history==1
                % Diffusion state history during addtime
                ds_add_h=[[0,S_DS];[dt_addtime,sm{diff_state_idx}]];
            end

        else % In that case we have to subdivide
            n_substeps=ceil(dt_addtime/max_dt);
            dt_subaddtime=dt_addtime/n_substeps;
            sub_xy=zeros(n_substeps+1,2); % define the intermediate positions sub_x and sub_y
            sub_xy(1,:)=[x_h, y_h]; % starting position
            c_diff_state=S_DS; % Current diffusion state

            if record_ds_history==1
                ds_add_h=nan(n_substeps+1,2); % Diffusion state history during addtime
                ds_add_h(1,:)=[0,S_DS]; % Assign first state to starting state at time 0
            end

            for k=1:n_substeps
                % Look for a potential change in diffusion state
                [sm{n_diff_state_idx}, sm{n_sp_idx}]=get_pattern_transition([sub_xy(k,1),sub_xy(k,2)], ...
                    V_ras(c_diff_state), D_ras(c_diff_state), c_diff_state,K,C,dt_subaddtime,sm_par,im_par); % Get new diffusion coefficient

                % Update position and check if a potential change in diffusion state is realized or not
                [sub_xy(k+1,1),sub_xy(k+1,2),sm]=get_new_XY_pct(sub_xy(k,1),sub_xy(k,2), sm, ...
                    D_ras, dt_subaddtime, S, sm_par, im_par);

                % Reassign current diffusion state
                c_diff_state=sm{diff_state_idx};

                if record_ds_history==1
                    % Diffusion state history during addtime
                    ds_add_h(k+1,:)=[k*dt_subaddtime,c_diff_state];
                end
            end
            x_h_d_tmp=sub_xy(end,1); % ending position
            y_h_d_tmp=sub_xy(end,2);
        end

        N_DS=sm{diff_state_idx}; % Reassign starting diffusion state

    else % Nothing needs to be changed
        x_h_d_tmp=x_h;
        y_h_d_tmp=y_h;
        N_DS=S_DS; % Do no change the starting diffusion state

        if record_ds_history==1
            % Diffusion state history during addtime
            ds_add_h=[];
        end
    end

    %Then evaluate how it moves during frametime

    % First do a quick evaluation to estimate the diffused distance
    % Look for a potential change in diffusion state
    [sm{n_diff_state_idx}, sm{n_sp_idx}]=get_pattern_transition([x_h_d_tmp,y_h_d_tmp], ...
        V_ras(N_DS), D_ras(N_DS), N_DS,K,C,dt_frametime,sm_par,im_par);

    % Update position and check if a potential change in diffusion state is realized or not
    [x_h_d_tmp2, y_h_d_tmp2,sm_tmp2]=get_new_XY_pct(x_h_d_tmp, y_h_d_tmp, ...
        sm, D_ras, dt_frametime, S, sm_par, im_par);

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

        for k=1:n_substeps
            % Look for a potential change in diffusion state
            [sm{n_diff_state_idx}, sm{n_sp_idx}]=get_pattern_transition([sub_xy(k,1),sub_xy(k,2)], ...
                V_ras(c_diff_state), D_ras(c_diff_state), c_diff_state,K,C,dt_subframetime,sm_par,im_par); % Get new diffusion coefficient

            % Update position and check if a potential change in diffusion state is realized or not
            [sub_xy(k+1,1),sub_xy(k+1,2),sm]=get_new_XY_pct(sub_xy(k,1),sub_xy(k,2), sm, ...
                D_ras, dt_subframetime, S, sm_par, im_par);

            % Reassign current diffusion state
            c_diff_state=sm{diff_state_idx};

            if record_ds_history==1
                % Diffusion state history during frametime
                ds_fr_h(k+1,:)=[dt_addtime+k*dt_subframetime,c_diff_state];
            end
        end
        x_h_d=sub_xy(end,1); % ending position
        y_h_d=sub_xy(end,2);
    else % if moved distance was small enough, or if there is no issue with fast exchange rates keep it as is
        x_h_d=x_h_d_tmp2; % Keep initial diffused position
        y_h_d=y_h_d_tmp2;
        sm=sm_tmp2; % Keep initial updated sm

        if record_ds_history==1
            % Diffusion state history during frametime
            ds_fr_h=[[dt_addtime,N_DS];[dt_addtime+dt_frametime,sm{diff_state_idx}]];
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

%% case of no diffuse PSF
else 

    % Look for a change in diffusion state
    [sm{n_diff_state_idx}, sm{n_sp_idx}]=get_pattern_transition([x_h,y_h], V_ras(S_DS), D_ras(S_DS), S_DS,K,C,dt,sm_par,im_par); % Get new diffusion coefficient

    % Update position and check if a potential change in diffusion state is realized or not
    [x_h_d, y_h_d,sm]=get_new_XY_pct(x_h, y_h, sm, D_ras, dt, S, sm_par, im_par);

    if record_ds_history==1
        ds_h=[[0,S_DS];[dt_addtime+dt_frametime,sm{diff_state_idx}]]; % Full diffusion state history during addtime + frametime
    end
end

if record_ds_history==1
    sm{diff_state_trace_idx}=vertcat(sm{diff_state_trace_idx}, ds_h);
else
    sm{diff_state_trace_idx}=horzcat(sm{diff_state_trace_idx}, [current_frame;sm{diff_state_idx}]);
end

% Update coordinates on high resolution image
sm{x_idx}=x_h_d;
sm{y_idx}=y_h_d;
% Update sub_coordinates on high resolution image if
% diffuse_psf
if use_diffuse_psf==1
    if diffused_distance>diffuse_psf_radius
        sm{sub_x_idx}=sub_xy(:,1);
        sm{sub_y_idx}=sub_xy(:,2);
    else
        sm{sub_x_idx}=[];
        sm{sub_y_idx}=[];
    end
end

% Update diffusion track for next frame
sm{x_track_idx}=vertcat(sm{x_track_idx}, [x_h_d,current_frame+1]);
sm{y_track_idx}=vertcat(sm{y_track_idx}, [y_h_d,current_frame+1]);


