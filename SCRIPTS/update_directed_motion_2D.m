function [p_h_out, v_h_out, p_h_out_np]=update_directed_motion_2D(p_h, v_h, dt, c_ds, n_ds, sm_par, sm_pattern_indices, im_par)

% PURPOSE:
% Update_directed_motion in 2D for a single mol with initial speed V
% The main idea is to update the velocity vector first based on the persistence 
% length of the current pattern, and then check whether this new velocity vector is 
% compatible with the current speed of the molecule. If, with this new direction 
% and speed the molecule ends up out of the pattern or out of the image, then 
% another approach needs to be taken. We then take for the direction the one 
% which is closest to the initial direction, and that allows for the molecule to 
% stay within the pattern at its current speed. In the case the molecule goes out
% of the image, we simply reverse the velocity vector and the molecule is supposed 
% to bounce back. If a molecule was positioned on an isolated piece of pattern 
% that doesn't allow movement with the current speed, then the molecule is assigned a speed of zero.
%
% INPUTS:
% 	p_h: the xy-coordinate on high-resolution image [raster]
%   v_h: the xy-coordinate of the initial velocity vector on high-resolution image [raster]
%	dt: the evolution time in [seconds]
%   c_ds: current sub_pattern id (or diffusion state) of the molecule
%   n_ds: potential new sub_pattern id (or diffusion state) of the molecule
%	sm_par: the sm parameters
%   sm_pattern_indices: indices of virtual sample subpatterns
%	im_par: the imaging parameters
%
% OUTPUTS:
%	p_h_out: [pixels] the updated xy-coordinate on high-resolution image
%	v_h_out: [pixels.s-1] the updated velocity vector on high-resolution image
%	p_h_out_np: [pixels] the updated xy-coordinate on high-resolution image if molecule to jump to new pattern

% MODIFICATION HISTORY:
%	D.Bourgeois, November 2020: version > simulate_palm_vsn16.3
%	D.Bourgeois, January 2021:  Update for change of diffusion state within single pattern
%	D.Bourgeois, February 2023, introduce sm_pattern_indices, now disconnected from sm_par

cp_id=sm_par.n_sp_id==sm_par.D_confined(c_ds); % Pattern id for current diffusion state
np_id=sm_par.n_sp_id==sm_par.D_confined(n_ds); % Pattern id for new diffusion state
c_w_pattern=sm_pattern_indices.w_patterns(cp_id).w;

show = 0; % Set to 1 for debug; 2 for final view

%% Define the orientation of the velocity vector based on the persistence length of the current pattern
% Get the coordinates of the init-dir circle around current point
x_r=sm_par.V_init_dir(c_ds).x;
y_r=sm_par.V_init_dir(c_ds).y;

% place on main image
x=p_h(1)+x_r;
y=p_h(2)+y_r;

% make sure circle is inside image
sz=im_par.binning*[im_par.n, im_par.m];
x(x<1)=1+sm_par.margin_factor;
x(x>sz(1))=sz(1)-sm_par.margin_factor;
y(y<1)=1+sm_par.margin_factor;
y(y>sz(2))=sz(2)-sm_par.margin_factor;

%look at the intersection with the current pattern
ind = sub2ind(sz,round(x),round(y)); % the index associated to the circle
ind_cp=intersect(ind,c_w_pattern); % the indices of the circle on the current pattern

if show==1 % Set up image if we want to see it
    I=zeros(sz);
    I(c_w_pattern)=10;
    I(round(p_h(1)), round(p_h(2)))=50;
    I(ind)=20;
    r2=1e+3*sm_par.persistence_length(c_ds)/im_par.raster*im_par.binning;
    D=2*max([1,ceil(r2)])+1; % The search diameter has to exceed what can be accessed in reality
else
    I=[];
    D=[];
end

%First treat case where the sm can reach a point within pattern
%at a distance corresponding to the persistence length
if ~isempty(ind_cp)  
    % Get the possible directions (use round values of x y z
    % position to avoid nonsymmetric angular values over full circle)
    [u_m, I]=get_directions_on_2D_pattern(ind_cp, round(p_h(1)), round(p_h(2)), sz, show, I, D);
    
    %Treat case of fully isotropic choice: in that case keep current direction
    if norm(u_m(:,1))<1e-3 % Do not use a limit of zero, the dispersion might not be exactly 0 for a perfect circle
        if size(u_m,2)~=1
            error('Isotropic direction with more than 1 velocity group is not possible !');
        end
        v_h_out=v_h;
    else % General case of an anisotropic distribution
        if size(u_m,2)==1 % This situation might happen when a molecule bump into a new pattern. In that case, we will choose 
            % a random direction as compared to the current one           
            ind_f=ind_cp(max([1,round(rand*numel(ind_cp))])); % Pick randomly one of the values
            [x_f, y_f]=ind2sub(sz,ind_f); % x,y coordinates corresponding to tp_f
            
            % Define the chosen unit vector u_f
            v=[x_f'-round(p_h(1)); y_f'-round(p_h(2))];
            u_f=v./vecnorm(v); 
        else % In that case we have several directions possible. Typically the case of a molecule moving along a microtubule or a 1D filament
            %The norm of the velocity vectors are
            u_m_n=vecnorm(u_m);
            %The angles between current and new possible velocity vectors are
            dth=acos(dot(repmat(v_h/norm(v_h),[size(u_m,2),1])',u_m./u_m_n));
            
            if any(isnan(dth)) % Deal with molecules that are stuck at the border of the image, 
                %That is remove the vectors that have zero length and end up with nan values for dth
                u_m=u_m(:,~isnan(dth));
                dth=dth(~isnan(dth));
            end
            
            %Let's define the probability law for dth, as an exponentially decaying
            %function of the square of dth over the range pi/4 (empirical)
            if ~isempty(dth)
                th_p=(u_m_n.^sm_par.dispersion_selectivity).*exp(-abs(dth.^2)/(pi/4)^2);
                th_p=th_p/sum(th_p); % Normalize to 1
                w_u=find(cumsum(th_p)>=rand,1); % The chosen direction is the first one for which the cumulative sum exceeds a random number between 0 an 1
                u_f=u_m(:,w_u)/u_m_n(w_u); % Normalize the vector
            else % This should never happen as in this case the molecule should automatically turn around, but just in case
                disp('Molecule may have hit image border ! ');
                warning('Setting speed to opposite current speed for current molecule !');
                v_h_out=-v_h; % Simply set the speed to zero
                u_f=[];
            end
        end
        
        %Get the proposed velocity vector
        if ~isempty(u_f)
            v_h_out=u_f'*1e+3*sm_par.V(c_ds)/im_par.raster*im_par.binning; % in [raster.s-1]
        end
    end
    
    %Update the position 
    p_h_out=p_h+dt*v_h_out;

    %Now we need to check that the new position is within the pattern.
    %Indeed, the direction was calculated based on the persistence length of the pattern, that is a 
    %possibly large radius from sm_par.V_init_dir,
    %but the actual radius is from sm_par.V_circle and depending on local curvature the proposed velocity 
    %vector might end up out of the pattern. and if not we choose the closest position available
    
    %First check that the molecule ends up within the field of view
    if ~(any([p_h_out(1)<1 p_h_out(2)<1]) || any([p_h_out(1)>sz(1) p_h_out(2)>sz(2)])) % check if molecule is inside the FOV
        
        check_ind = sub2ind(sz,round(p_h_out(1)),round(p_h_out(2))); % the index associated to the new position
        check_inside = ismember(check_ind,c_w_pattern);
        
        if check_inside==0 % In that case the molecule does not end up with in the pattern and 
            %we need to change the reached position !    
            %Select velocity vector which is closest to the initial direction, and that allows for the molecule to 
            % stay within the pattern at its current speed.
            [p_h_out, v_h_out]=get_p_v_2D_out(p_h, v_h, sm_par.V_circle(c_ds), sm_par.V(c_ds), dt, im_par.raster, im_par.binning, sm_par.margin_factor, c_w_pattern, sz);
        end
    else % In that case the molecule has got out of the field of view
        disp('Molecule may have hit image border ! ');
        warning('Setting speed to opposite current speed for current molecule !');
        v_h_out=-v_h; % Simply reverse the speed
        p_h_out=p_h; % Simply do not move the molecule
    end
    
else % Treat case where no point in pattern can be reached at specified velocity
    %In that case, we take the same approach as above, and will look for the least deviation from current velocity vector
    [p_h_out, v_h_out]=get_p_v_2D_out(p_h, v_h, sm_par.V_circle(c_ds), sm_par.V(c_ds), dt, im_par.raster, im_par.binning, sm_par.margin_factor, c_w_pattern, sz);
end

%Draw the velocity on I
if show==1
    figure(1)
    clf
    line([p_h(2),p_h_out(2)],[p_h(1),p_h_out(1)],'Color','red','LineWidth',3);
    hold on
    imagesc(I)
    axis image
    hold off
    drawnow
    set(gcf,'color','white')
    title('Chosen direction')
    input('Ok ?','s')
end



%Finally check if the current velocity vector crosses a possible new pattern if a transition is allowed.
% if n_ds~=c_ds
if n_ds~=c_ds  && any(cp_id~=np_id) % Corrected 17/01/21
    n_w_pattern=sm_pattern_indices.w_patterns(np_id).w;
    
    % In that case get the intersection of the current velocity vector with the new pattern
    n_pos=round(norm(dt*v_h));
    on_new_pattern=false;
    k=n_pos;
   
    while ~on_new_pattern && (k>0)
        p_h_out_np=p_h+k/n_pos*dt*v_h; % Possible position on the new pattern
        ind_np = sub2ind(sz,round(p_h_out_np(1)),round(p_h_out_np(2))); % the index associated to the new position
        on_new_pattern = ismember(ind_np,n_w_pattern);
        k=k-1;
    end
    %Only consider p_h_out_np if the new pattern can be reached !
    if on_new_pattern~=1
        p_h_out_np=[];
    end
else
    p_h_out_np=[];
end
