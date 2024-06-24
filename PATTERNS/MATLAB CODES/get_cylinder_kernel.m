grid_size=256;
cylinder_length=64; % Size of the ellipsoid radius on X axis [in pixels)
cylinder_b=50; % Size of the ellipsoid radius on Y axis [in pixels)
% cylinder_c=200; % Size of the ellipsoid radius on Z axis [in pixels)
cylinder_c=50; % Size of the ellipsoid radius on Z axis [in pixels)
cylinder_thick=1; % Thickness of the ellipsoid [in pixels) (This will be approximate)
z_limit=110; %50; % Overall height of object [in pixels]
theta=0; % angle between cylinder axis and x,y planeclf

show=1;
show_surface=1;

save_matfile=0;
matfile='D:\Mes Documents\MATLAB\SIMULATION\PALM\SOFTWARE\PALM_SIMULATIONS_vsn14.4\PATTERNS\Cylinder_complete_Z.mat';

[x,y,z]=meshgrid(1:grid_size,1:grid_size,1:grid_size);
x=x-grid_size/2;
y=y-grid_size/2;
z=z-grid_size/2;

% Here is the equation of the 3D volume
inner_radius=1-cylinder_thick/cylinder_b;
outer_radius=1+cylinder_thick/cylinder_b;
t=x>(-cylinder_length/2) & x<(cylinder_length/2) & (sqrt((y./cylinder_b).^2+(z./cylinder_c).^2) < outer_radius) & (sqrt((y./cylinder_b).^2+(z./cylinder_c).^2) > inner_radius) & z>(-z_limit/2) & z<(z_limit/2);

t2=permute(t,[1,3,2]); % permute y and z
t2=imrotate(t2,theta,'bilinear','crop'); % rotation around z
t2=permute(t2,[1,3,2]); % permute y and z
if show==1
    figure(1)
    clf
    h=vol3d('cdata', t2*100);
    alphamap('default')
    colormap('jet')
    axis image
    view(3);
    set(gcf,'color','white')
end

if save_matfile==1
    kernel_pattern=t2;
    save(matfile,'kernel_pattern');
end

if show_surface==1
    figure(2)
    colormap gray
    step=10;
    n_steps=fix(grid_size/step);
    ij_steps=fix(sqrt(n_steps));
    k=1;
    for i=1:ij_steps
        for j=1:ij_steps
            
            subplot(ij_steps,ij_steps,ij_steps*(i-1)+j)
            k=k+step;
            imagesc(~t2(:,:,k)); axis image
            title(['Z=',num2str(k)]);
        end
    end
end

%Z projection
figure(3)
colormap gray
t3=sum(t2,3);
imagesc(t3); axis image

