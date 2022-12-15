
close all
grid_size_xy=256;
grid_size_z=64;

%Bacterial pattern
bead_thick=10; % Thickness raster size over which beads might be positionned 
bead_height=5; % Height at which beads might be positionned 

show=1;
show_z0_section=1;
save_matfile=0;
matfile='C:\Users\bourgeoi\Documents\MATLAB\SIMULATION\PALM\SOFTWARE\PALM_SIMULATIONS_vsn16.3\PATTERNS\Calibration_Beads_TIRF.mat';

% End input
im=zeros(grid_size_xy,grid_size_xy,grid_size_z);
im(:,:,round(bead_height+1-bead_thick/2):round(bead_height+1+bead_thick/2)-1)=1;

if show==1
    figure(1)
    clf
    h=vol3d('cdata', im);
    alphamap('default')
    colormap('jet')
    axis image
    view(3);
end

if show_z0_section==1
    figure(2)
    im0=im(:,:,round(grid_size_z/2));
    imagesc(im0)
    axis image
    colormap('gray')
end


if save_matfile==1
    kernel_pattern=im;
    save(matfile,'kernel_pattern');
end

