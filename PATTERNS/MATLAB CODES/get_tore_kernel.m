grid_size=64;
tore_inner_radius=8; % Size of the tore
tore_height=2; % Thickness of the tore. Total diameter will be 2*(tore_inner_radius+tore_height)
theta=0; % angle between tore axis and x,y plane

show=1;
save_matfile=0;
matfile='D:\Mes Documents\MATLAB\SIMULATION\PALM\PALM_SIMULATIONS_vsn14\PATTERNS\Flat_Z_ring.mat';

[x,y,z]=meshgrid(1:grid_size,1:grid_size,1:grid_size);
x=x-grid_size/2;
y=y-grid_size/2;
z=z-grid_size/2;
t=((tore_inner_radius-sqrt(x.^2+y.^2)).^2+z.^2)<tore_height^2;
t2=permute(t,[1,3,2]); % permute y and z
t2=imrotate(t2,theta,'bilinear','crop'); % rotation around z
t2=permute(t2,[1,3,2]); % permute y and z
if show==1
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

