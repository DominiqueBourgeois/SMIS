function [im,Proj]=get_SMIS_3D_Microtubules(par)


% show=1;

x_dim=par.x_dim;
y_dim=par.y_dim;
z_dim=par.z_dim;

N=par.N; % # of microtubules

mt_length=par.mt_length; % Microtubule length [in pixels)
mt_radius=par.mt_radius; % Microtubule radius [in pixels)
mt_thick=par.mt_thick; % Microtubule thickness [in pixels)
mt_rot=par.mt_rot; % Microtubule rotation [in degrees)

% mt_length=500; 
% mt_radius=2; % Size of the cylinder radius on Y axis [in pixels)
% mt_thick=1; % Thickness of the cylinder[in pixels) (This will be approximate) Use mt_thick < 2*mt_radius the attempt getting an hollow cylinder

cell_border_offset=par.cell_border_offset; % [Pixels]

% Shift relative to optical axis
move_sample_to_coverslip=par.move_sample_to_coverslip;
z_offset=par.z_offset;

id_mt=1; % Pattern id associated to mt;

save_projection_image=1;

%%

grid_size=[x_dim, y_dim, z_dim];

%Shift
shift_val_XY= grid_size(1)/2;
shift_val_Z=0;
shift_X=round(shift_val_XY*(rand(N,1)-0.5)); % angle between cylinder axis and x,y planeclf
shift_Y=round(shift_val_XY*(rand(N,1)-0.5)); % angle between cylinder axis and x,y planeclf
shift_Z=round(shift_val_Z*(rand(N,1)-0.5)); % angle between cylinder axis and x,y planeclf

%Rotation
rot_val_XY=90; % [Deg]
rot_val_Z=mt_rot; % [Deg]
theta_XY=rot_val_XY*(rand(N,1)-0.5); % angle between cylinder axis and x,y planeclf
theta_XZ=rot_val_Z*(rand(N,1)-0.5); % angle between cylinder axis and x,y planeclf
theta_YZ=rot_val_Z*(rand(N,1)-0.5); % angle between cylinder axis and x,y planeclf

[x0,y0,z0]=meshgrid(1:grid_size(1),1:grid_size(2),1:grid_size(3));
im = x0<0;

x0=x0-grid_size(1)/2;
y0=y0-grid_size(2)/2;
z0=z0-grid_size(3)/2;
inner_radius=1-mt_thick/mt_radius;
outer_radius=1+mt_thick/mt_radius;

MyWaitBar = waitbar(0,'Generating microtubules ...');


for k=1:N
    disp(['Computing MT #: ', num2str(k)]);
    waitbar(k/N,MyWaitBar);

    % Here is the equation of the 3D volume
    x=x0-shift_X(k);
    y=y0-shift_Y(k);
    z=z0-shift_Z(k);
    
    t=x>(-mt_length/2) & x<(mt_length/2) & (sqrt((y./mt_radius).^2+(z./mt_radius).^2) < outer_radius) ...
        & (sqrt((y./mt_radius).^2+(z./mt_radius).^2) > inner_radius) & z>(-z_dim/2) & z<(z_dim/2);
    
    %XY rotation
    t2=permute(t,[2,1,3]); % permute y and z
    t2=imrotate(t2,theta_XY(k),'bilinear','crop'); % rotation around z
    t2=permute(t2,[2,1,3]); % permute y and z
    
    %XZ rotation
    t2=permute(t2,[3,2,1]); % permute y and z
    t2=imrotate(t2,theta_XZ(k),'bilinear','crop'); % rotation around z
    t2=permute(t2,[3,2,1]); % permute y and z
    
    %YZ rotation
    t2=permute(t2,[1,3,2]); % permute y and z
    t2=imrotate(t2,theta_YZ(k),'bilinear','crop'); % rotation around z
    t2=permute(t2,[1,3,2]); % permute y and z
       
    im=t2 | im; % Merge the kernel
end
close(MyWaitBar)

im=im*id_mt;

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

title('Virtual 3D Microtubules')

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


