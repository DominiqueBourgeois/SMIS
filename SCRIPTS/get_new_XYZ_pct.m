function [x_h_d, y_h_d, z_h_d, sm] = get_new_XYZ_pct(x_h, y_h, z_h, sm, D_ras, dt, S, sm_par, im_par)

% PURPOSE:
% Get new position for a molecule initially at x,y diffusing in 2D for time dt with diffusion coeff D
%
% INPUTS:
%	x_h: the x-coordinate on high-resolution image
%	y_h: the y-coordinate on high-resolution image
%	z_h: the z-coordinate on high-resolution image
%	sm: the current sm
%   D_ras: diffusion coefficients in raster units [raster]^2/s-1, both in XY and in Z
%	dt: the evolution time in [seconds]
%   S: size of a high-resolution image
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%	x_h_d: the moved x-coordinate on high-resolution image
%	y_h_d: the moved y-coordinate on high-resolution image
%	z_h_d: the moved z-coordinate on high-resolution image
%   sm: updated sm
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15
%	D.Bourgeois, November 2019: version > simulate_palm_vsn15.3
%	D.Bourgeois, December 2019: version > simulate_palm_vsn15.3
%	D.Bourgeois, September 2020: version > simulate_palm_vsn16.2 (prevent_loosing_molecules)
%	D.Bourgeois, November 2020: version > simulate_palm_vsn16.3 (directed diffusion)
%	D.Bourgeois, January 2021: Added possibility to change of diffusion state within single pattern
%	D.Bourgeois, April 2021: Added possibility have state's change independent of diffusion
%	D.Bourgeois, September 2022, optimized for parallel computing

%Extract the useful indices
% x_idx=1;
% y_idx=2;
% z_idx=3;
% x_track_idx=4;
% y_track_idx=5;
% z_track_idx=6;
% sub_x_idx=7;
% sub_y_idx=8;
% sub_z_idx=9;
v_x_idx=10;
v_y_idx=11;
v_z_idx=12;
c_sp_idx=13;
n_sp_idx=14;
% id_idx=15;
% bleached_idx=16;
% activated_idx=17;
diff_state_idx=18;
n_diff_state_idx=19;
% diff_state_trace_idx=20;
% matched_idx=21;

w=sm_par.w_patterns; %	w: indices of patterns in high resolution image

c_sp=sm{c_sp_idx}; % current subpattern
c_ds=sm{diff_state_idx}; % current diffusion state
n_sp=sm{n_sp_idx}; % potential new subpattern
n_ds=sm{n_diff_state_idx}; % potential new current diffusion state

D_ras_XY=D_ras(1,:); % The diffusion coeffs in XY (raster units)
D_ras_Z=D_ras(2,:); % The diffusion coeffs in Z (raster units)

%Handle Diff independant transitions
if sm_par.DIT(c_ds)==1
    force_change=1;
else
    force_change=0;
end

if force_change==1
    %Search for the closest pixel in new pattern and move molecule
    %there
    if n_sp~=c_sp
        P = move_molecule_to_new_pattern([x_h,y_h,z_h],w(sm_par.n_sp_id==n_sp).w,im_par);
        x_h_d=P(1)+rand;
        y_h_d=P(2)+rand;
        z_h_d=P(3)+rand;
        sm{diff_state_idx}=n_ds; % update diffusion state
        sm{c_sp_idx}=n_sp; % update current pattern
        if sm_par.V(n_ds)>0 % % update velocity
            v_h = get_initial_directions(x_h_d, y_h_d, z_h_d, n_sp, n_ds, im_par, sm_par);
            sm{v_x_idx}=v_h(1);
            sm{v_y_idx}=v_h(2);
            sm{v_z_idx}=v_h(3);
        end
        return % Stop here !
    end
end

%Handle directed motion if necessary
if sm_par.V(c_ds)>0
    [p_h_out,v_h_out,p_h_out_np]=update_directed_motion_3D_pct([x_h,y_h,z_h],[sm{v_x_idx}, sm{v_y_idx}, sm{v_z_idx}], dt, c_ds, n_ds, sm_par,im_par);
    if isempty(p_h_out_np) % Normal case, no change of pattern possible
        x_h=p_h_out(1);
        y_h=p_h_out(2);
        z_h=p_h_out(3);
        sm{v_x_idx}=v_h_out(1);
        sm{v_y_idx}=v_h_out(2);
        sm{v_z_idx}=v_h_out(3);
        
    else % Here a change of pattern can be envisaged
        % Consider the potential new position in the new pattern
        x_h=p_h_out_np(1);
        y_h=p_h_out_np(2);
        z_h=p_h_out_np(3);
        % And reset the velocity
        sm{v_x_idx}=[];
        sm{v_y_idx}=[];
        sm{v_z_idx}=[];
        
        %Apply diffusion: If the single molecule stays in the new pattern,
        %keep this solution.
        x_h_d=x_h+sqrt(2*D_ras_XY(c_ds)*dt)*randn;
        y_h_d=y_h+sqrt(2*D_ras_XY(c_ds)*dt)*randn;
        z_h_d=z_h+sqrt(2*D_ras_Z(c_ds)*dt)*randn;

        %Prevent molecule to diffuse too far away by 10% of image size
        if im_par.prevent_diffusion_out_FOV==1 % 1 if diffusion out of more than 10% of FOV size is prevented==1
            if x_h_d<(-0.1*S(1))
                x_h_d=-0.1*S(1);
            end
            if y_h_d<(-0.1*S(2))
                y_h_d=-0.1*S(2);
            end
            if z_h_d<(-0.1*S(3))
                z_h_d=-0.1*S(3);
            end
            if x_h_d>(1.1*S(1))
                x_h_d=1.1*S(1);
            end
            if y_h_d>(1.1*S(2))
                y_h_d=1.1*S(2);
            end
            if z_h_d>(1.1*S(3))
                z_h_d=1.1*S(3);
            end
        end
        
        %Check that the molecule stays within the field of view

        if ~(any([x_h_d<1 y_h_d<1 z_h_d<1]) || any([x_h_d>S(1) y_h_d>S(2) z_h_d>S(3)])) % check if molecule is inside the FOV
            w_b = sub2ind(S,round(x_h_d),round(y_h_d),round(z_h_d)); % the associated index
            if min(abs(w_b-w(sm_par.n_sp_id==sm_par.D_confined(n_ds)).w))==0 % Success, molecule jumped to new pattern !
                jump_ok=1;
            else
                jump_ok=0;
            end
        elseif n_sp==0 % If molecule is out of the Field-of-view but Is floating in the background
            % In this case This is ok, let the molecule diffuse away
            jump_ok=1;
        else % This case is when there is rapid diffusion and the molecule is on a non-background pattern
            jump_ok=0;
        end
        if jump_ok==1
            sm{c_sp_idx}=n_sp; % update current pattern
            sm{diff_state_idx}=n_ds; % update current diffusion state
            % If there is directed motion, the velocity must be reset and reinitialized when there is a change of pattern
            v_h = get_initial_directions(x_h_d, y_h_d, z_h_d, n_sp, n_ds, im_par, sm_par);
            sm{v_x_idx}=v_h(1);
            sm{v_y_idx}=v_h(2);
            sm{v_z_idx}=v_h(3);
            return
        else % The jump to a new pattern was unsuccessful, keep going on the current pattern !
            x_h=p_h_out(1);
            y_h=p_h_out(2);
            z_h=p_h_out(3);
            sm{v_x_idx}=v_h_out(1);
            sm{v_y_idx}=v_h_out(2);
            sm{v_z_idx}=v_h_out(3);
        end
        %If the single molecule does not stay in the new pattern after applying diffusion, go back to the original scheme
    end
end

%Treat special case of D_ras=0 + avoids getting trapped
if D_ras_XY(c_ds)==0 && D_ras_Z(c_ds)==0
    x_h_d=x_h; % Do not change x and y positions
    y_h_d=y_h;
    z_h_d=z_h;
    
    if n_ds~=c_ds && sm_par.D_confined(n_ds)==sm_par.D_confined(c_ds) % Case where there is
        % a change in diffusion state within a single subpattern
        sm{diff_state_idx}=n_ds; % eventually update current diffusion state
        if sm_par.V(c_ds)>0 % eventually update current velocity Changed 17/01/21
            v_h=[sm{v_x_idx}, sm{v_y_idx}, sm{v_z_idx}]; % This is the current velocity vector
            % Simply scale the velocity vector to the new speed
            sm{v_x_idx}=v_h(1)/norm(v_h)*1e+3*sm_par.V(n_ds)/im_par.raster*im_par.binning;
            sm{v_y_idx}=v_h(2)/norm(v_h)*1e+3*sm_par.V(n_ds)/im_par.raster*im_par.binning;
            sm{v_z_idx}=v_h(3)/norm(v_h)*1e+3*sm_par.V(n_ds)/im_par.raster*im_par.binning;
        else % In that case there was no specific direction for speed, we need to initialize it.
            v_h = get_initial_directions(x_h_d, y_h_d, z_h_d, n_sp, n_ds, im_par, sm_par);
            sm{v_x_idx}=v_h(1);
            sm{v_y_idx}=v_h(2);
            sm{v_z_idx}=v_h(3);
        end
    end
    return
end

%Now handle diffusion
on_pattern=false;

while ~on_pattern
    x_h_d=x_h+sqrt(2*D_ras_XY(c_ds)*dt)*randn;
    y_h_d=y_h+sqrt(2*D_ras_XY(c_ds)*dt)*randn;
    z_h_d=z_h+sqrt(2*D_ras_Z(c_ds)*dt)*randn;
   
    %Prevent molecule to diffuse too far away by 10% of image size
    if im_par.prevent_diffusion_out_FOV==1 % 1 if diffusion out of more than 10% of FOV size is prevented==1
        if x_h_d<(-0.1*S(1))
            x_h_d=-0.1*S(1);
        end
        if y_h_d<(-0.1*S(2))
            y_h_d=-0.1*S(2);
        end
        if z_h_d<(-0.1*S(3))
            z_h_d=-0.1*S(3);
        end
        
        if x_h_d>(1.1*S(1))
            x_h_d=1.1*S(1);
        end
        if y_h_d>(1.1*S(2))
            y_h_d=1.1*S(2);
        end
        if z_h_d>(1.1*S(3))
            z_h_d=1.1*S(3);
        end
    end
    
    if ~(any([x_h_d<1 y_h_d<1 z_h_d<1]) || any([x_h_d>S(1) y_h_d>S(2) z_h_d>S(3)])) % check if molecule is inside the FOV
        w_b = sub2ind(S,round(x_h_d),round(y_h_d),round(z_h_d)); % the associated index
        
        if ismember(w_b,w(sm_par.n_sp_id==sm_par.D_confined(c_ds)).w) % the indices of current pattern is that where sm_par.n_sp_id=sm_par.D_confined(c_ds) CHANGE 16 01 2021
            on_pattern=true;  % molecule stays on the current pattern, we're ok.
            if n_ds~=c_ds && sm_par.D_confined(n_ds)==sm_par.D_confined(c_ds) % Case where there is
                % a change in diffusion state within a single subpattern
                sm{diff_state_idx}=n_ds; % update current diffusion state
                if sm_par.V(c_ds)>0 % eventually update current velocity Changed 17/01/21
                    v_h=[sm{v_x_idx}, sm{v_y_idx}, sm{v_z_idx}]; % This is the current velocity vector
                    % Simply scale the velocity vector to the new speed
                    sm{v_x_idx}=v_h(1)/norm(v_h)*1e+3*sm_par.V(n_ds)/im_par.raster*im_par.binning;
                    sm{v_y_idx}=v_h(2)/norm(v_h)*1e+3*sm_par.V(n_ds)/im_par.raster*im_par.binning;
                    sm{v_z_idx}=v_h(3)/norm(v_h)*1e+3*sm_par.V(n_ds)/im_par.raster*im_par.binning;
                else % In that case there was no specific direction for speed, we need to initialize it.
                    v_h = get_initial_directions(x_h_d, y_h_d, z_h_d, n_sp, n_ds, im_par, sm_par);
                    sm{v_x_idx}=v_h(1);
                    sm{v_y_idx}=v_h(2);
                    sm{v_z_idx}=v_h(3);
                end
            end
        elseif n_sp~=c_sp % look if there is a possible change in subpattern.
            if min(abs(w_b-w(sm_par.n_sp_id==sm_par.D_confined(n_ds)).w))==0
                on_pattern=true;  % molecule moved to the new pattern, we're ok.
                sm{c_sp_idx}=n_sp; % update current pattern
                sm{diff_state_idx}=n_ds; % update current diffusion state
                % now update current velocity vector if new diffusion state has an associated speed
                % the velocity must be reset and reinitialized when there is a change of pattern
                % if sm_par.V(c_ds)>0
                if sm_par.V(n_ds)>0 % Changed 17/01/21
                    v_h = get_initial_directions(x_h_d, y_h_d, z_h_d, n_sp, n_ds, im_par, sm_par);
                    sm{v_x_idx}=v_h(1);
                    sm{v_y_idx}=v_h(2);
                    sm{v_z_idx}=v_h(3);
                end
            end
        end
    else % If molecule is out of the Field-of-view
        if c_sp==0 % In this case This is ok, let the molecule diffuse away
            on_pattern=true;  % molecule is out of FOV but considered to be on the pattern, we're ok.
        end
    end
    
end




