function [Sample,Proj]=get_SMIS_3DCubesSpheres(par)


n=par.x_dim; %number of column
m=par.y_dim; %number of slices
nz=par.z_dim; %number of slices

is_cube=par.is_cube;

qPALM_option=par.qPALM;

%spacing between lines [pixel]
s_xy=par.s_xy; % size or diameter in xy
s_z=par.s_z; % size or diameter in z

n_row=par.n_row; % # of rows along x
n_col=par.n_col; % # of columns along y
n_col_z=par.n_col_z; % # of columns along z

r_xy_shift=par.random_xy_shift; % Random shift of positions
r_z_shift=par.random_z_shift; % Random shift of positions

%safety border
cell_border_offset=par.cell_border_offset; % [Pixels]

feature_id=1;

% Shift relative to optical axis
move_sample_to_coverslip=par.move_sample_to_coverslip;
z_offset=par.z_offset;

save_projection_image=1;

%%
im = double(zeros(m,n,nz));

if (n_col>n_row && m<n) || (n_col<n_row && m>n)
    tmp=n_row;
    n_row=n_col;
    n_col=tmp;
end

inc_r=fix(n/n_row)+1; % # of increments along rows
inc_c=fix(m/n_col)+1; % # of increments along columns
inc_c_z=fix(nz/n_col_z)+1; % # of increments along height

x0=inc_r:inc_r:inc_r*n_row;
y0=inc_c:inc_c:inc_c*n_col;
z0=inc_c_z:inc_c_z:inc_c_z*n_col_z;

offset_x=(x0(1)+x0(end)-(n-1))/2;
offset_y=(y0(1)+y0(end)-(m-1))/2;
offset_z=(z0(1)+z0(end)-(nz-1))/2;

x0=x0-fix(offset_x);
y0=y0-fix(offset_y);
z0=z0-fix(offset_z);

% x0=x0-round(inc_r/2);
% y0=y0-round(inc_c/2);
% z0=z0-round(inc_c_z/2);

v=1; % pattern id, to be increased for qPALM

[x,y,z]=meshgrid(1:n,1:m,1:nz);
disp('Creating pattern ...')
for i=1:n_row
    x_c=round(x0(i));
    for j=1:n_col
        y_c=round(y0(j));
        for k=1:n_col_z
            z_c=round(z0(k));

            if r_xy_shift>0
                x_c=x_c+r_xy_shift*(rand-0.5);
                y_c=y_c+r_xy_shift*(rand-0.5);
            end
            if r_z_shift>0
                z_c=z_c+r_z_shift*(rand-0.5);
            end

%             if is_cube==1
%                 if mod(s_xy,2)==1 && mod(s_z,2)==1
%                     t=abs(x-x_c)<=fix(s_xy/2) & abs(y-y_c)<=fix(s_xy/2) & abs(z-z_c)<=fix(s_z/2);
%                 elseif mod(s_xy,2)==1 && mod(s_z,2)~=1
%                     t=abs(x-x_c)<=fix(s_xy/2) & abs(y-y_c)<=fix(s_xy/2) & z<(z_c+s_z/2+1) & z>=(z_c-s_z/2+1);
%                 elseif mod(s_xy,2)~=1 && mod(s_z,2)==1
%                     t=x<(x_c+s_xy/2+1) & x>=(x_c-s_xy/2+1) & y<(y_c+s_xy/2+1) & y>=(y_c-s_xy/2+1) & abs(z-z_c)<=fix(s_z/2);
%                 else
%                     t=x<(x_c+s_xy/2+1) & x>=(x_c-s_xy/2+1) & y<(y_c+s_xy/2+1) & y>=(y_c-s_xy/2+1) & z<(z_c+s_z/2+1) & z>=(z_c-s_z/2+1);
%                 end
%             else
%                 t=((x-x_c).^2/s_xy^2+(y-y_c).^2/s_xy^2+(z-z_c).^2/s_z^2)<=1;
%             end

            if is_cube==1
                if mod(s_xy,2)==1 && mod(s_z,2)==1
                    t=abs(x-x_c)<=fix(s_xy/2) & abs(y-y_c)<=fix(s_xy/2) & abs(z-z_c)<=fix(s_z/2);
                elseif mod(s_xy,2)==1 && mod(s_z,2)~=1
                    t=abs(x-x_c)<=fix(s_xy/2) & abs(y-y_c)<=fix(s_xy/2) & z<(z_c+s_z/2 ) & z>=(z_c-s_z/2 );
                elseif mod(s_xy,2)~=1 && mod(s_z,2)==1
                    t=x<(x_c+s_xy/2 ) & x>=(x_c-s_xy/2 ) & y<(y_c+s_xy/2 ) & y>=(y_c-s_xy/2 ) & abs(z-z_c)<=fix(s_z/2);
                else
                    t=x<(x_c+s_xy/2 ) & x>=(x_c-s_xy/2 ) & y<(y_c+s_xy/2 ) & y>=(y_c-s_xy/2 ) & z<(z_c+s_z/2 ) & z>=(z_c-s_z/2 );
                end
            else
                t=((x-x_c).^2/s_xy^2+(y-y_c).^2/s_xy^2+(z-z_c).^2/s_z^2)<=1;
            end

            if qPALM_option==1
                im(t==1)=v;
                v=v+1;
            else
                im=im | t;
            end
        end
    end
end
disp('Done !')

im=im*feature_id;

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


%%
disp('Done !');


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

if is_cube==1
    title('Virtual 3D Cubes')
else
    title('Virtual 3D Spheres')
end

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

