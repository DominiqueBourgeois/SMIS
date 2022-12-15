% function t = get_tore_kernel(grid_size, tore_inner_radius, tore_height, theta, show)

grid_size=256;
ellipsoid_a=64; % Size of the ellipsoid radius on X axis [in pixels)
ellipsoid_b=24; % Size of the ellipsoid radius on Y axis [in pixels)
ellipsoid_c=24; % Size of the ellipsoid radius on Z axis [in pixels)
ellipsoid_thick=2; % Thickness of the ellipsoid [in pixels)
theta=0; % angle between ellipsoid axis and x,y plane

show=1;
save_matfile=0;
matfile='D:\Mes Documents\MATLAB\SIMULATION\PALM\SOFTWARE\PALM_SIMULATIONS_vsn14.4\PATTERNS\Ellipsoid_test.mat';

[x,y,z]=meshgrid(1:grid_size,1:grid_size,1:grid_size);
x=x-grid_size/2;
y=y-grid_size/2;
z=z-grid_size/2;

% Here is the equation of the 3D volume
inner_radius=1-ellipsoid_thick/2;
outer_radius=1+ellipsoid_thick/2;
t=(sqrt((x./ellipsoid_a).^2+(y./ellipsoid_b).^2+(z./ellipsoid_c).^2) < outer_radius) & (sqrt((x./ellipsoid_a).^2+(y./ellipsoid_b).^2+(z./ellipsoid_c).^2) > inner_radius);

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
end

if save_matfile==1
    kernel_pattern=t2;
    save(matfile,'kernel_pattern');
end

