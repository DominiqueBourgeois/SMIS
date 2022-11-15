function [Sample,Proj]=get_SMIS_3D_EColi(par)

%Pattern of elliptical bacterium with periplasm

% show=1;

% cytoplasm_a=80; % Size of the cytoplasm radius on X axis [in pixels)
% cytoplasm_b=20; % Size of the cytoplasm radius on Y axis [in pixels)

cytoplasm_a=par.cytoplasm_a; % Size of the cytoplasm radius on X axis [in pixels)
cytoplasm_b=par.cytoplasm_b; % Size of the cytoplasm radius on Y axis [in pixels)
cytoplasm_c=par.cytoplasm_c; % Size of the cytoplasm radius on Z axis [in pixels)

include_nucleoid=par.include_nucleoid; % 1 if nucleoid to be included
nucleoid_a=par.nucleoid_a; % Size of the nucleoid radius on X axis [in pixels)
nucleoid_b=par.nucleoid_b; % Size of the nucleoid radius on Y axis [in pixels)
nucleoid_c=par.nucleoid_c; % Size of the nucleoid radius on Z axis [in pixels)

% periplasm_thick=2; % Thickness of the periplasm [in pixels)
periplasm_thick=par.periplasm_thick; % Thickness of the periplasm [in pixels)

n_cells_x=par.n_cells_y; % # of cells X dim
n_cells_y=par.n_cells_x; % # of cells Y dim

% cell_border_offset=100; % [Pixels]
% cell_spacing=5; % [Pixels]

cell_border_offset=par.cell_border_offset; % [Pixels]
cell_spacing=par.cell_spacing; % [Pixels]

% rotate_cells=1; % Set to 1 to randomly rotate cells ** CAREFUL THIS CAN PRODUCE WEIRD PIXEL VALUES ***
rotate_cells=par.rotate_cells; % Set to 1 to randomly rotate cells ** CAREFUL THIS CAN PRODUCE WEIRD PIXEL VALUES ***
xy_rot=20;
z_rot=180;

% Shift relative to optical axis
move_sample_to_coverslip=par.move_sample_to_coverslip;
z_offset=par.z_offset;
z_dim=par.z_dim;

id_outside=0; % Pattern id associated to elsewhere;
id_cytoplasm=1; % Pattern id associated to cytoplasm;
id_periplasm=2; % Pattern id associated to cytoplasm;
id_nucleoid=3; % Pattern id associated to nucleoid;

save_projection_image=1;

%%
clf
colormap('jet');

grid_size_xy=cell_spacing+max([cytoplasm_a,cytoplasm_b,cytoplasm_c])+periplasm_thick;
grid_size_z=z_dim;

%Set parameters
a0=cytoplasm_a-cytoplasm_b;
r0b=cytoplasm_b/2;
r0c=cytoplasm_c/2;


[x,y,z]=meshgrid(1:grid_size_xy,1:grid_size_xy,1:grid_size_z);
x=x-grid_size_xy/2;
y=y-grid_size_xy/2;
z=z-grid_size_z/2;

full_im=zeros(grid_size_xy*n_cells_x,grid_size_xy*n_cells_y,grid_size_z);
%% Create the object

x0=-a0/2;

% Draw cytoplasm
cytoplasm_center=(((y./r0b).^2+(z./r0c).^2 <= 1) & (x >= -a0/2 & x <= a0/2));

cytoplasm_left=((((x-x0)./r0b).^2+(y./r0b).^2+(z./r0c).^2 <= 1) & x<=x0);

x0=a0/2;
cytoplasm_right=((((x-x0)./r0b).^2+(y./r0b).^2+(z./r0c).^2 <= 1) & x>=x0);

cytoplasm=cytoplasm_center | cytoplasm_left | cytoplasm_right;
% cytoplasm=cytoplasm_center;

% Draw nucleoid
if include_nucleoid==1
    nucleoid=((x./(nucleoid_a/2)).^2+(y./(nucleoid_b/2)).^2+(z./(nucleoid_c/2)).^2 <= 1);
end

% Draw periplasm
inner_radius_2=1;
outer_radius_2=1+periplasm_thick/mean(0.5*[cytoplasm_b,cytoplasm_c]);


periplasm_center=(((y./r0b).^2+(z./r0c).^2 > inner_radius_2^2) & ((y./r0b).^2+(z./r0c).^2 <= outer_radius_2^2) & (x >= -a0/2 & x <= a0/2));

x0=-a0/2;
periplasm_left=((((x-x0)./r0b).^2+(y./r0b).^2+(z./r0c).^2 > inner_radius_2^2) & ...
    (((x-x0)./r0b).^2+(y./r0b).^2+(z./r0c).^2 <= outer_radius_2^2) & x<=x0);


x0=a0/2;
periplasm_right=((((x-x0)./r0b).^2+(y./r0b).^2+(z./r0c).^2 > inner_radius_2^2) & ...
    (((x-x0)./r0b).^2+(y./r0b).^2+(z./r0c).^2 <= outer_radius_2^2) & x>=x0);


periplasm=periplasm_center | periplasm_left | periplasm_right;

% Assemble
if include_nucleoid==1
    im=cytoplasm*id_cytoplasm+periplasm*id_periplasm+nucleoid*(id_nucleoid-id_cytoplasm);
    im(im~=id_cytoplasm & im~=id_periplasm & im~=id_nucleoid)=id_outside;
else
    im=cytoplasm*id_cytoplasm+periplasm*id_periplasm;
    im(im~=id_cytoplasm & im~=id_periplasm)=id_outside;
end

s=size(im);

if rotate_cells==1
    im_start=im;
    cell_rotation_x=xy_rot*rand(1,n_cells_x*n_cells_y); % [degrees] angle of the cell ** CAREFUL THIS CAN PRODUCE WEIRD PIXEL VALUES ***
    cell_rotation_y=xy_rot*rand(1,n_cells_x*n_cells_y); % [degrees] angle of the cell ** CAREFUL THIS CAN PRODUCE WEIRD PIXEL VALUES ***
    cell_rotation_z=z_rot*rand(1,n_cells_x*n_cells_y); % [degrees] angle of the cell ** CAREFUL THIS CAN PRODUCE WEIRD PIXEL VALUES ***
end

%% Create all cells
MyWaitBar = waitbar(0,'Generating cells ...');

for i=1:n_cells_x
    for j=1:n_cells_y
        waitbar(((i-1)*n_cells_y+j)/(n_cells_y*n_cells_x),MyWaitBar);

        clc
        disp(['Generating cell # : ', num2str(j+(i-1)*n_cells_y)])

        if rotate_cells==1

            rx=cell_rotation_x(j+n_cells_y*(i-1));
            ry=cell_rotation_z(j+n_cells_y*(i-1));
            rz=cell_rotation_y(j+n_cells_y*(i-1));

            im=im_start;

            if rx~=0
                tmp_im=permute(im,[1,3,2]); % permute y and z
                tmp_im=imrotate(tmp_im,rx,'bilinear','crop'); % rotation around z
                im=permute(tmp_im,[1,3,2]); % permute y and z
            end

            if ry~=0
                tmp_im=permute(im,[2,1,3]); % permute y and z
                tmp_im=imrotate(tmp_im,ry,'bilinear','crop'); % rotation around z
                im=permute(tmp_im,[2,1,3]); % permute y and z
            end

            if rz~=0
                tmp_im=permute(im,[3,2,1]); % permute y and z
                tmp_im=imrotate(tmp_im,rz,'bilinear','crop'); % rotation around z
                im=permute(tmp_im,[3,2,1]); % permute y and z
            end

            %clean up
            im=round(im); % avoid interpolated values


            %Clean up borders of bacteria: a cytoplasm pixel should not
            %contact directly the outside
            tmp_im=im;
            if periplasm_thick>0
                w=find(im==id_cytoplasm);
                for k=1:numel(w)

                    [x_,y_,z_] = ind2sub(s,w(k)); % the associated index

                    if (z_-1)<1 || (z_+1)>size(im,3) || (y_-1)<1 || (y_+1)>size(im,2)
                        warndlg('Image is too small, increase Z-dimension !')
                        Sample=[];
                        Proj=[];
                        close(MyWaitBar)
                        return
                    else
                        si1=im(x_-1,y_-1:y_+1,z_-1:z_+1);
                        si2=im(x_+1,y_-1:y_+1,z_-1:z_+1);
                        si3=im(x_,y_-1,z_-1:z_+1);
                        si4=im(x_,y_+1,z_-1:z_+1);
                        si5=im(x_,y_,z_-1:2:z_+1);

                        around=[si1(:);si2(:);si3(:);si4(:);si5(:)];

                        if any(around==id_outside)
                            tmp_im(w(k))=id_outside;
                        end
                    end
                end
            end
            im=tmp_im;

            %Clean up borders of bacteria: a periplasm pixel should not
            %contact directly the nucleoid
            if include_nucleoid==1
                tmp_im=im;
                w=find(im==id_periplasm);
                for k=1:numel(w)

                    [x_,y_,z_] = ind2sub(s,w(k)); % the associated index

                    if (z_-1)<1 || (z_+1)>size(im,3) || (y_-1)<1 || (y_+1)>size(im,2)
                        warndlg('Image is too small, increase Z-dimension !')
                        Sample=[];
                        Proj=[];
                        close(MyWaitBar)
                        return
                    else
                        si1=im(x_-1,y_-1:y_+1,z_-1:z_+1);
                        si2=im(x_+1,y_-1:y_+1,z_-1:z_+1);
                        si3=im(x_,y_-1,z_-1:z_+1);
                        si4=im(x_,y_+1,z_-1:z_+1);
                        si5=im(x_,y_,z_-1:2:z_+1);

                        around=[si1(:);si2(:);si3(:);si4(:);si5(:)];

                        if any(around==id_nucleoid)
                            tmp_im(w(k))=id_cytoplasm;
                        end
                    end
                end
                im=tmp_im;
            end


        end

        full_im((i-1)*grid_size_xy+1:i*grid_size_xy,(j-1)*grid_size_xy+1:j*grid_size_xy,:)=im;
    end
end
close(MyWaitBar)


%Move along optical axis
if move_sample_to_coverslip==1
    w=find(full_im>0);
    [~,~,z_]=ind2sub(size(full_im),w);
    z_offset=min(z_);
    full_im=circshift(full_im,-z_offset+1,3);
elseif z_offset<0
    full_im(:,:,1:-z_offset)=0;
    full_im=circshift(full_im,z_offset,3);
elseif z_offset>0
    full_im(:,:,end-z_offset:end)=0;
    full_im=circshift(full_im,z_offset,3);
end

Sample=zeros(size(full_im,1)+2*cell_border_offset,size(full_im,2)+2*cell_border_offset,size(full_im,3));
Sample(cell_border_offset+1:cell_border_offset+size(full_im,1),...
    cell_border_offset+1:cell_border_offset+size(full_im,2),:)=full_im;


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

title('Virtual 3D E.coli Cells')

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


