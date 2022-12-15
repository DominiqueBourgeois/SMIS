function v_h = get_initial_directions(x_h, y_h, z_h, sp, ds, im_par, sm_par)

% PURPOSE:
% Get initial velocity vectors for sm's at x,y,z diffusing in 2D or 3D at
% speed V
%
% INPUTS:
%	x_h: the x-coordinates on high-resolution image [raster]
%	y_h: the y-coordinates on high-resolution image [raster]
%	z_h: the z-coordinates on high-resolution image [raster] in case of 3D
%   sp: the current subpattern ids for the sm's
%   ds: the current diffusion states for the sm's
%	im_par: the imaging parameters
%	sm_par: the sm parameters
%
% OUTPUTS:
%	v_h: the initial speed vectors on high-resolution image [raster.s-1]
%
% MODIFICATION HISTORY:
%	D.Bourgeois, November 2020: version > simulate_palm_vsn16.3

w_patterns=sm_par.w_patterns;
show = 0; % Set to 1 for debug; 2 for final view

%2D case
if im_par.simul_3D==0
    %Go through all molecules
    v_h=zeros(numel(x_h),2); % Define the velocity vectors for all sm
    for i=1:numel(x_h)
        c_sp=sp(i); % current subpattern
        c_sp_id=find(sm_par.n_sp_id==c_sp,1);
        c_ds=ds(i); % current diffusion state % Changed 17/01/21
        %         if sm_par.V(c_sp_id)>0
        if sm_par.V(c_ds)>0 % Changed 17/01/21
            % radius to look around [pixels] in high res image, according
            % to set persistence length
            %             x=sm_par.V_init_dir(c_sp_id).x;
            %             y=sm_par.V_init_dir(c_sp_id).y;
            x=sm_par.V_init_dir(c_ds).x; % Changed 17/01/21
            y=sm_par.V_init_dir(c_ds).y;
            
            % place on main image
            x=x_h(i)+x;
            y=y_h(i)+y;
            
            % make sure disk is inside image
            sz=im_par.binning*[im_par.n, im_par.m];
            x(x<1)=1+sm_par.margin_factor;
            x(x>sz(1))=sz(1)-sm_par.margin_factor;
            y(y<1)=1+sm_par.margin_factor;
            y(y>sz(2))=sz(2)-sm_par.margin_factor;
            
            %look at the intersection
            ind = sub2ind(sz,round(x),round(y)); % the index associated to the sphere
            ind_cp=intersect(ind,w_patterns(c_sp_id).w); % the indices of the sphere on the current pattern
            
            if show==1 % Set up image if we want to see it
                I=zeros(sz);
                I(w_patterns(c_sp_id).w)=10;
                I(round(x_h(i)), round(y_h(i)))=50;
                I(ind)=20;
            else
                I=[];
            end
            
            if ~isempty(ind_cp)
                %The next two lines are only useful if show=1
                %                 r2=1e+3*sm_par.persistence_length(c_sp_id)/im_par.raster*im_par.binning;
                r2=1e+3*sm_par.persistence_length(c_ds)/im_par.raster*im_par.binning; % Changed 17/01/21
                D=2*max([1,ceil(r2)])+1; % The search diameter has to exceed what can be accessed in reality
                
                [u_mean, I]=get_directions_on_2D_pattern(ind_cp, x_h(i), y_h(i), sz, show, I, D);
                
                if show==1
                    colormap jet
                    xc=round(x_h(i));
                    yc=round(y_h(i));
                    I_loc=I(max([1,xc-D]):min([xc+D,sz(1)]),max([1,yc-D]):min([yc+D,sz(2)]));
                    figure(2)
                    clf
                    imagesc(I_loc);
                    axis image
                    set(gcf,'color','white')
                    title('2D Local view of the pattern plus search ball plus selected regions')
                end
                
                %Now calculate the final chosen angle
                
                %Treat case of fully isotropic choice: in that case choose a
                %random direction
                if norm(u_mean)==0
                    if size(u_mean,2)~=1
                        error('Isotropic direction with more than 1 velocity group is not possible !');
                    end
                    u_f=rand(2,1); u_f=u_f/norm(u_f); % Choose a random orientation factor
                else % General case
                    if size(u_mean,2)==1
                        u_f=u_mean/norm(u_mean);
                    else
                        %The probabilities to choose a particular direction is set to be proportional
                        %to the norm of the unit vectors(large dispersion = large channel)
                        th_p=vecnorm(u_mean).^sm_par.dispersion_selectivity;
                        th_p=th_p/sum(th_p); % Normalize to 1
                        u_f=u_mean(:,find(cumsum(th_p)>=rand,1)); % The chosen direction is the first one for which the cumulative sum exceeds a random number between 0 an 1
                        u_f=u_f/norm(u_f); % Normalize to one
                    end
                end
                
                %Get the final velocity vector for the sm #i
                %                 v_h(i,:)=u_f'*1e+3*sm_par.V(c_sp_id)/im_par.raster*im_par.binning; % in [raster.s-1]
                v_h(i,:)=u_f'*1e+3*sm_par.V(c_ds)/im_par.raster*im_par.binning; % in [raster.s-1] % Changed 17/01/21
                
                %Draw the velocity on I
                if show==1
                    figure(3)
                    clf
                    imagesc(I_loc);
                    axis image
                    set(gcf,'color','white')
                    title('Chosen direction')
                    hold on
                    line([size(I_loc,1)/2,size(I_loc,1)/2+u_f(2)*D/2],[size(I_loc,2)/2,size(I_loc,2)/2+ u_f(1)*D/2],'Color','red','LineWidth',3);
                    hold off
                    drawnow
                    input('Ok ?','s')
                end
            else % Treat case where no point in pattern can be reached at specified velocity
                v_h(i,:)=0; % Simply set the speed to zero
            end
        else % No velocity to define
            v_h(i,:)=0;
        end
        
    end
    
    if show==2 % Get a final global view
        figure(1)
        clf
        I=zeros(sz);
        I(w_patterns(c_sp_id).w)=10;
        hold on
        u_h=v_h'./vecnorm(v_h');
        for i=1:numel(x_h)
            I(round(x_h(i)), round(y_h(i)))=50;
            line([y_h(i),y_h(i)+u_h(2,i)*D/2],[x_h(i),x_h(i)+ u_h(1,i)*D/2],'Color','red','LineWidth',3);
        end
        imagesc(I);
        axis image
        hold off
        drawnow
        set(gcf,'color','white')
        title('Chosen direction')
        input('Ok ?','s')
    end
elseif im_par.simul_3D==1 % 3D Case
    %Go through all molecules
    v_h=zeros(numel(x_h),3); % Define the velocity vectors for all sm
    for i=1:numel(x_h)
        c_sp=sp(i); % current subpattern
        c_sp_id=find(sm_par.n_sp_id==c_sp,1);
        c_ds=ds(i); % current diffusion state % Changed 17/01/21
        %         if sm_par.V(c_sp_id)>0
        if sm_par.V(c_ds)>0 % Changed 17/01/21
            % radius to look around [pixels] in high res image, according
            % to set persistence length
            %             x=sm_par.V_init_dir(c_sp_id).x;
            %             y=sm_par.V_init_dir(c_sp_id).y;
            %             z=sm_par.V_init_dir(c_sp_id).z;
            x=sm_par.V_init_dir(c_ds).x; % Changed 17/01/21
            y=sm_par.V_init_dir(c_ds).y;
            z=sm_par.V_init_dir(c_ds).z;
            
            % place on main image
            x=x_h(i)+x;
            y=y_h(i)+y;
            z=z_h(i)+z;
            
            % make sure disk is inside image
            sz=im_par.binning*[im_par.n, im_par.m, im_par.nz];
            x(x<1)=1+sm_par.margin_factor;
            x(x>sz(1))=sz(1)-sm_par.margin_factor;
            y(y<1)=1+sm_par.margin_factor;
            y(y>sz(2))=sz(2)-sm_par.margin_factor;
            z(z<1)=1+sm_par.margin_factor;
            z(z>sz(3))=sz(3)-sm_par.margin_factor;
            
            %look at the intersection
            ind = sub2ind(sz,round(x),round(y),round(z)); % the index associated to the sphere
            ind_cp=intersect(ind,w_patterns(c_sp_id).w); % the indices of the sphere on the current pattern
            
            if show==1 % Set up image if we want to see it
                I=zeros(sz);
                I(w_patterns(c_sp_id).w)=10;
                I(round(x_h(i)), round(y_h(i)), round(z_h(i)))=50;
                I(ind)=20;
            else
                I=[];
            end
            
            if ~isempty(ind_cp)
                %The next two lines are only useful if show=1
                %                 r2=1e+3*sm_par.persistence_length(c_sp_id)/im_par.raster*im_par.binning;
                r2=1e+3*sm_par.persistence_length(c_ds)/im_par.raster*im_par.binning; % Changed 17/01/21
                D=2*max([1,ceil(r2)])+1; % The search diameter has to exceed what can be accessed in reality
                
                [u_mean, I]=get_directions_on_3D_pattern(ind_cp, x_h(i), y_h(i), z_h(i), sz, show, I, D);
                
                if show==1
                    colormap jet
                    xc=round(x_h(i));
                    yc=round(y_h(i));
                    zc=round(z_h(i));
                    I_loc=I(max([1,xc-D]):min([xc+D,sz(1)]),max([1,yc-D]):min([yc+D,sz(2)]),max([1,zc-D]):min([zc+D,sz(3)]));
                    figure(2)
                    clf
                    vol3d('cdata', I_loc);
                    alphamap('default')
                    axis image
                    view(3);
                    set(gcf,'color','white')
                    title('3D Local view of the pattern plus search ball plus selected regions')
                    
                    figure(3)
                    step=1;
                    n_steps=fix(size(I_loc,3)/step);
                    ij_steps=ceil(sqrt(n_steps));
                    loc_k=0;
                    for loc_i=1:ij_steps
                        for loc_j=1:ij_steps
                            loc_k=loc_k+step;
                            if loc_k<=size(I_loc,3)
                                subplot(ij_steps,ij_steps,ij_steps*(loc_i-1)+loc_j)
                                imagesc(I_loc(:,:,loc_k)); axis image
                                caxis([0, max(I_loc(:))]);
                                title(['Z=',num2str(loc_k)]);
                            end
                        end
                    end
                    sgtitle('Cross-sections through image in figure 2')
                end
                
                %Now calculate the final chosen angle
                
                %Treat case of fully isotropic choice: in that case choose a
                %random direction
                if norm(u_mean)==0
                    if size(u_mean,2)~=1
                        error('Isotropic direction with more than 1 velocity group is not possible !');
                    end
                    u_f=rand(3,1); u_f=u_f/norm(u_f); % Choose a random orientation factor
                else % General case
                    if size(u_mean,2)==1
                        u_f=u_mean/norm(u_mean);
                    else
                        %The probabilities to choose a particular direction is set to be proportional
                        %to the norm of the unit vectors(large dispersion = large channel)
                        th_p=vecnorm(u_mean).^sm_par.dispersion_selectivity;
                        th_p=th_p/sum(th_p); % Normalize to 1
                        u_f=u_mean(:,find(cumsum(th_p)>=rand,1)); % The chosen direction is the first one for which the cumulative sum exceeds a random number between 0 an 1
                        u_f=u_f/norm(u_f); % Normalize to one
                    end
                end
                
                %Get the final velocity vector for the sm #i
                %                 v_h(i,:)=u_f'*1e+3*sm_par.V(c_sp_id)/im_par.raster*im_par.binning; % in [raster.s-1]
                v_h(i,:)=u_f'*1e+3*sm_par.V(c_ds)/im_par.raster*im_par.binning; % in [raster.s-1] % Changed 17/01/21
                
                %Draw the velocity on I
                if show==1
                    figure(4)
                    clf
                    vol3d('cdata', I_loc);
                    alphamap('default')
                    axis image
                    view(3);
                    set(gcf,'color','white')
                    title('Chosen direction')
                    hold on
                    line([size(I_loc,1)/2,size(I_loc,1)/2+u_f(2)*D/2],[size(I_loc,2)/2,size(I_loc,2)/2+ u_f(1)*D/2], [size(I_loc,3)/2,size(I_loc,3)/2+ u_f(3)*D/2],'Color','red','LineWidth',3);
                    hold off
                    drawnow
                    input('Ok ?','s')
                end
            else % Treat case where no point in pattern can be reached at specified velocity
                v_h(i,:)=0; % Simply set the speed to zero
            end
        else % No velocity to define
            v_h(i,:)=0;
        end
        
    end
    
    if show==2 % Get a final global view
        figure(1)
        clf
        I=zeros(sz);
        I(w_patterns(c_sp_id).w)=10;
        hold on
        u_h=v_h'./vecnorm(v_h');
        for i=1:numel(x_h)
            I(round(x_h(i)), round(y_h(i)), round(z_h(i)))=50;
            line([y_h(i),y_h(i)+u_h(2,i)*D/2],[x_h(i),x_h(i)+ u_h(1,i)*D/2], [z_h(i),z_h(i)+ u_h(3,i)*D/2],'Color','red','LineWidth',3);
        end
        vol3d('cdata', I);
        alphamap('default')
        axis image
        view(3);
        hold off
        drawnow
        set(gcf,'color','white')
        title('Chosen direction')
        input('Ok ?','s')
    end
end




