function [im,Proj]=get_SMIS_3D_NPCs(par)


% show=1;

x_dim=par.y_dim;
y_dim=par.x_dim;
z_dim=par.z_dim;

N=par.N; % # of NPCs
qPALM_option=par.qPALM; %Set to 1 to assign a different pixel value to all nup96s in the field for qPALM experiments

% nup96_size_xy=1;  % [pixels] Nup96 size (area is l^2)
% nup96_size_z=1;  % [pixels] Nup96 size (area is l^2)
% npc_in=5; % [pixels] NPC inner radius
% npc_out=7; % [pixels] NPC outer radius
% npc_thickness=6; % [pixels] NPC height

nup96_size_xy=par.nup96_size_xy;  % [pixels] Nup96 size (area is l^2)
nup96_size_z=par.nup96_size_z;  % [pixels] Nup96 size (area is l^2)
npc_in=par.npc_in; % [pixels] NPC inner radius
npc_out=par.npc_out; % [pixels] NPC outer radius
npc_thick=par.npc_thick; % [pixels] NPC height


% min_sep=10; % [pixels] Minimum separation between two NPC's
min_sep=par.min_sep; % [pixels] Minimum separation between two NPC's

allow_rotation=par.allow_rotation; % Set to 1 if rotation allowed
rot_xy=par.rot_xy;  % XY rotation
rot_z=par.rot_z;  % Z rotation

cell_border_offset=par.cell_border_offset; % [Pixels]

% Shift relative to optical axis
move_sample_to_coverslip=par.move_sample_to_coverslip;
z_offset=par.z_offset;

id_npc=1; % Pattern id associated to mt;

save_projection_image=1;

%%

npc_rotation_x=rot_xy*(rand(1,N)-0.5); % [degrees] angle of the cell ** CAREFUL THIS CAN PRODUCE WEIRD PIXEL VALUES ***
npc_rotation_z=rot_z*(rand(1,N)-0.5); % [degrees] angle of the cell ** CAREFUL THIS CAN PRODUCE WEIRD PIXEL VALUES ***
npc_rotation_y=rot_xy*(rand(1,N)-0.5); % [degrees] angle of the cell ** CAREFUL THIS CAN PRODUCE WEIRD PIXEL VALUES ***

% npc_im_size=3*npc_out+1;
npc_im_size=2*npc_out+1;


% Choose the closest 2^n number
loc_n=1;
while npc_im_size>2^loc_n
    loc_n=loc_n+1;
end
npc_im_size=2^loc_n;

npc_im=zeros(npc_im_size,npc_im_size,npc_im_size);

lower_layer_position=round(npc_im_size/2-npc_thick/2);
upper_layer_position=round(npc_im_size/2+npc_thick/2);

lower_layer=ceil(lower_layer_position-nup96_size_z/2):fix(lower_layer_position+nup96_size_z/2);
upper_layer=ceil(upper_layer_position-nup96_size_z/2):fix(upper_layer_position+nup96_size_z/2);

% lower_layer=round(lower_layer_position-nup96_size_z/2):round(lower_layer_position+nup96_size_z/2);
% upper_layer=round(upper_layer_position-nup96_size_z/2):round(upper_layer_position+nup96_size_z/2);


v=1;
%Inner ring
for theta=0:2*pi/8:(2*pi-2*pi/8)
    x=npc_im_size/2+npc_in*cos(theta); x=round(x);
    y=npc_im_size/2+npc_in*sin(theta); y=round(y);
    
    npc_im(x:x+nup96_size_xy-1, y:y+nup96_size_xy-1,lower_layer)=v;
    if qPALM_option==1
        v=v+1;
    end
    npc_im(x:x+nup96_size_xy-1, y:y+nup96_size_xy-1,upper_layer)=v;
    if qPALM_option==1
        v=v+1;
    end
    
end
%Outer ring
for theta=0:2*pi/8:(2*pi-2*pi/8)
    x=npc_im_size/2+npc_out*cos(theta); x=round(x);
    y=npc_im_size/2+npc_out*sin(theta); y=round(y);
    npc_im(x:x+nup96_size_xy-1, y:y+nup96_size_xy-1,lower_layer)=v;
    if qPALM_option==1
        v=v+1;
    end
    npc_im(x:x+nup96_size_xy-1, y:y+nup96_size_xy-1,upper_layer)=v;
    if qPALM_option==1
        v=v+1;
    end   
end


%% Full image
im3D=zeros(x_dim,y_dim,z_dim); % The full 3D image
im2D=zeros(x_dim,y_dim); % Mask to ensure that the NPC's are well separated

sep=round(min_sep+2*npc_out); % That's the minimum distance between two NPC centers
border=max([cell_border_offset,sep]); % To avoid problems at image borders

n2=x_dim-2*border;
m2=y_dim-2*border;

%Define stocheometry
if qPALM_option==1
    V=v-1;
else
    V=1;
end

n_pix=numel(find(npc_im>0 & npc_im<=V));


MyWaitBar = waitbar(0,'Generating NPCs ...');

% Place each individual NPC
for i=1:N
    waitbar(i/N,MyWaitBar);

    %Eventually rotate
    if allow_rotation==1
        rx=npc_rotation_x(i);
        ry=npc_rotation_z(i);
        rz=npc_rotation_y(i);
        
        %Rotation around x
        tmp_im=permute(npc_im,[1,3,2]); % permute y and z
        tmp_im=imrotate(tmp_im,rx,'bilinear','crop'); % rotation around z
        tmp_im=permute(tmp_im,[1,3,2]); % permute y and z
        
        %Rotation around y
        tmp_im=permute(tmp_im,[2,1,3]); % permute y and z
        tmp_im=imrotate(tmp_im,ry,'bilinear','crop'); % rotation around z
        tmp_im=permute(tmp_im,[2,1,3]); % permute y and z
        
        %Rotation around z
        tmp_im=permute(tmp_im,[3,2,1]); % permute y and z
        tmp_im=imrotate(tmp_im,rz,'bilinear','crop'); % rotation around z
        npc_im_rotated=permute(tmp_im,[3,2,1]); % permute y and z
        
        w=find(npc_im_rotated>0);
        [~,s]=sort(npc_im_rotated(w),'descend');
        npc_im_rotated(w(s(1:n_pix)))=1:V;
        npc_im_rotated(w(s(n_pix+1:end)))=0;
    else
        npc_im_rotated=npc_im;
    end
    
    
    sep_ok=0;
    trial_number=1;
    while ~sep_ok==1 && trial_number<1000
        X=randi(n2)+ border;
        Y=randi(m2)+ border;
        if im2D(X,Y)==0 && (X-npc_im_size/2)>0 && (X+npc_im_size/2-1)<=size(im3D,1) && (Y-npc_im_size/2)>0 && (Y+npc_im_size/2-1)<=size(im3D,2)
            sep_ok=1;
        end
        trial_number=trial_number+1;
    end

    if sep_ok==0
        warndlg('Could not place NPCs: increase image XY size !')
        im=[];
        Proj=[];
        return
    end

    im2D(X-sep:X+sep,Y-sep:Y+sep)=1; % Mask the corresponding 2D region

    z_min=(z_dim-size(npc_im_rotated,3))/2+1;
    z_max=z_min+size(npc_im_rotated,3)-1;

    if qPALM_option==0
        im3D(X-npc_im_size/2:X+npc_im_size/2-1,Y-npc_im_size/2:Y+npc_im_size/2-1,z_min:z_max)=npc_im_rotated; % Fill the corresponding 3D region
    else
        
        npc_im_rotated(npc_im_rotated>0)=npc_im_rotated(npc_im_rotated>0)+(i-1)*V;
        
        im3D(X-npc_im_size/2:X+npc_im_size/2-1,Y-npc_im_size/2:Y+npc_im_size/2-1,z_min:z_max)=npc_im_rotated; % Fill the corresponding 3D region
    end
end


close(MyWaitBar)

im=im3D*id_npc;

%Move along optical axis
if move_sample_to_coverslip==1
    w=find(im>0);
    [~,~,z_]=ind2sub(size(im),w);
    z_offset=min(z_);
    im=circshift(im,-z_offset+1,3);
elseif z_offset<0
    im(:,:,1:-z_offset)=0;
    im=circshift(im,z_offset,3);
elseif z_offset>0
    im(:,:,end-z_offset:end)=0;
    im=circshift(im,z_offset,3);
end

Sample=zeros(size(im,1)+2*cell_border_offset,size(im,2)+2*cell_border_offset,size(im,3));
Sample(cell_border_offset+1:cell_border_offset+size(im,1),...
    cell_border_offset+1:cell_border_offset+size(im,2),:)=im;


%% Show the cell
figure(1)
clf
colormap('jet');

disp('Displaying cells ...')
vol3d('cdata', Sample);
alphamap('default')
axis image
view(3);
xlabel('X [pixel]')
ylabel('Y [pixel]')
zlabel('Z [pixel]')

title('Virtual 3D NPCs')

if save_projection_image==1
    Proj=sum(Sample,3);
    figure(2)
    clf
    set(gcf,'Color','w')
    imagesc(Proj);
    xlabel('X [pixel]')
    ylabel('Y [pixel]')
    title('Z-projection')

    axis image
    colormap('gray')
else
    Proj=[];
end

disp('Done !');


