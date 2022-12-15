function [x,y,z] = get_sphere_coordinates(r)

% PURPOSE:
% Get x,y,z coordinates of sphere of radius r
%
% INPUTS:
%	r: the sphere radius [pixels]
%
% OUTPUTS:
%	x,y,z: the 2D cartesian coordinates [pixels]
%
% MODIFICATION HISTORY:
%	D.Bourgeois, November 2020: version > simulate_palm_vsn16.3

show=0;
show_surface=0;

D=2*max([1,ceil(r)])+1; % The search diameter has to exceed what can be accessed in reality

%define x,y values around x0, y0, z0 within radius r
[x,y,z]=meshgrid(1:D,1:D,1:D);

% recenter
x=x-(D+1)/2;
y=y-(D+1)/2;
z=z-(D+1)/2;

% keep the central circle
xyz_ok=(x.^2+y.^2+z.^2)<D^2/4 & (x.^2+y.^2+z.^2)>(D-2)^2/4; % Empirical formula: D-2 for the lower bounderies was found by trial and error

x=x(xyz_ok);
y=y(xyz_ok);
z=z(xyz_ok);

if show==1
    kernel_pattern=xyz_ok;
    figure(1)
    clf
    vol3d('cdata', kernel_pattern*100);
    alphamap('default')
    %     colormap('jet')
    axis image
    view(3);
    set(gcf,'color','white')
    
    
    if show_surface==1
        figure
        colormap gray
        step=1;
        n_steps=fix(D/step);
        ij_steps=ceil(sqrt(n_steps));
        k=0;
        for i=1:ij_steps
            for j=1:ij_steps
                k=k+step;
                if k<=size(kernel_pattern,3)
                    subplot(ij_steps,ij_steps,ij_steps*(i-1)+j)
                    imagesc(~kernel_pattern(:,:,k)); axis image
                end
                title(['Z=',num2str(k)]);
            end
        end
    end
end