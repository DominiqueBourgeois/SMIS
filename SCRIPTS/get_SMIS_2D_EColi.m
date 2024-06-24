function MySample=get_SMIS_2D_EColi(par)

%Pattern of elliptical bacterium with periplasm

% show=1;
% save_tiffile=0;
% tiffile='C:\Users\bourgeoi\Documents\MATLAB\SIMULATION\PALM\SOFTWARE\SMIS_vsn1.3\PATTERNS\2D\2D_RodShaped_Bacteria_Field.tif';

% cytoplasm_a=80; % Size of the cytoplasm radius on X axis [in pixels)
% cytoplasm_b=20; % Size of the cytoplasm radius on Y axis [in pixels)

cytoplasm_a=par.cytoplasm_a; % Size of the cytoplasm radius on X axis [in pixels)
cytoplasm_b=par.cytoplasm_b; % Size of the cytoplasm radius on Y axis [in pixels)

% nucleoid_a=25; % Size of the nucleoid radius on X axis [in pixels)
% nucleoid_b=5; % Size of the nucleoid radius on Y axis [in pixels)

include_nucleoid=par.include_nucleoid; % 1 if nucleoid to be included
nucleoid_a=par.nucleoid_a; % Size of the nucleoid radius on X axis [in pixels)
nucleoid_b=par.nucleoid_b; % Size of the nucleoid radius on Y axis [in pixels)


% periplasm_thick=2; % Thickness of the periplasm [in pixels)
periplasm_thick=par.periplasm_thick; % Thickness of the periplasm [in pixels)

% n_cells_x=10; % # of cells X dim
% n_cells_y=10; % # of cells Y dim

n_cells_x=par.n_cells_y; % # of cells X dim
n_cells_y=par.n_cells_x; % # of cells Y dim

% cell_border_offset=100; % [Pixels]
% cell_spacing=5; % [Pixels]

cell_border_offset=par.cell_border_offset; % [Pixels]
cell_spacing=par.cell_spacing; % [Pixels]

% rotate_cells=1; % Set to 1 to randomly rotate cells ** CAREFUL THIS CAN PRODUCE WEIRD PIXEL VALUES ***
rotate_cells=par.rotate_cells; % Set to 1 to randomly rotate cells ** CAREFUL THIS CAN PRODUCE WEIRD PIXEL VALUES ***

id_cytoplasm=1; % Pattern id associated to cytoplasm;
id_periplasm=2; % Pattern id associated to cytoplasm;
id_nucleoid=3; % Pattern id associated to nucleoid;
id_outside=0; % Pattern id associated to elsewhere;

% pixel_size=100; % [nm] SMIS Pixel size in nanometer
% binning=4; % SMIS binning

%%
clf
colormap('jet');

n_cells=n_cells_x*n_cells_y;

grid_size=cell_spacing+max([cytoplasm_a,cytoplasm_b])+periplasm_thick;


full_im=zeros(grid_size*n_cells_x,grid_size*n_cells_y);

if rotate_cells==1
    cell_rotation=90*rand(1,n_cells); % angle of the cell ** CAREFUL THIS CAN PRODUCE WEIRD PIXEL VALUES ***
else
    cell_rotation=0;
end


%Define 1 bacterium

[x,y]=meshgrid(1:grid_size,1:grid_size);
x=x-grid_size/2;
y=y-grid_size/2;

a0=cytoplasm_a-cytoplasm_b;
r0=cytoplasm_b/2;

% Draw cytoplasm
cytoplasm_c=(x >= -a0/2 & x <= a0/2 & y >= -r0 & y <= r0);

x0=-a0/2;
y0=0;
cytoplasm_l=((x-x0).^2+(y-y0).^2 <= r0^2 & x<=x0);

x0=a0/2;
y0=0;
cytoplasm_r=((x-x0).^2+(y-y0).^2 <= r0^2 & x>=x0);

cytoplasm=cytoplasm_c | cytoplasm_l | cytoplasm_r;

% Draw periplasm
r1=r0+periplasm_thick;

periplasm_c=((x >= -a0/2 & x <= a0/2) & ((y > r0 & y <= r1) |  (y < -r0 & y > -r1)));

x0=-a0/2;
y0=0;
periplasm_l=((x-x0).^2+(y-y0).^2 > r0^2 & (x-x0).^2+(y-y0).^2 <= r1^2 & x<=x0);

x0=a0/2;
y0=0;
periplasm_r=((x-x0).^2+(y-y0).^2 > r0^2 & (x-x0).^2+(y-y0).^2 <= r1^2 & x>=x0);

periplasm=periplasm_c | periplasm_l | periplasm_r;

% Draw nucleoid
if include_nucleoid==1
    nucleoid=(sqrt((x./nucleoid_a).^2+(y./nucleoid_b).^2) <= 1);

    im=cytoplasm*id_cytoplasm+periplasm*id_periplasm+nucleoid*(id_nucleoid-id_cytoplasm);
    im(im<id_cytoplasm & im<id_periplasm & im<id_nucleoid)=id_outside;
else
    im=cytoplasm*id_cytoplasm+periplasm*id_periplasm;
    im(im<id_cytoplasm & im<id_periplasm)=id_outside;
end

s=size(im);

for i=1:n_cells_x
    for j=1:n_cells_y
        clc
        disp(['Generating cell # : ', num2str(j+(i-1)*n_cells_y)])

        if rotate_cells==1
            single_bact=imrotate(im,cell_rotation(j+(i-1)*n_cells_y),'bilinear','crop'); % rotation around z
            single_bact=round(single_bact); % avoid interpolated values
            %Clean up border of bacteria: a cytoplasm pixel should not
            %contact directly the outside
            if periplasm_thick>0
                single_bact2=single_bact;
                w=find(single_bact==id_cytoplasm);
                for k=1:numel(w)
                    around=[single_bact(w(k)-1),single_bact(w(k)+1), ...
                        single_bact(w(k)-1+s(1)),single_bact(w(k)+s(1)),single_bact(w(k)+1+s(1)), ...
                        single_bact(w(k)-1-s(1)),single_bact(w(k)-s(1)),single_bact(w(k)+1-s(1))];
                    if any(around==id_outside)
                        single_bact2(w(k))=id_outside;
                    end
                end
                single_bact=single_bact2;
            end
                       
            
            %Clean up borders of bacteria: a periplasm pixel should not
            %contact directly the nucleoid
            if include_nucleoid==1
                single_bact2=single_bact;
                w=find(single_bact==id_periplasm);
                for k=1:numel(w)
                    around=[single_bact(w(k)-1),single_bact(w(k)+1), ...
                        single_bact(w(k)-1+s(1)),single_bact(w(k)+s(1)),single_bact(w(k)+1+s(1)), ...
                        single_bact(w(k)-1-s(1)),single_bact(w(k)-s(1)),single_bact(w(k)+1-s(1))];
                    if any(around==id_nucleoid)
                        single_bact2(w(k))=id_cytoplasm;
                    end
                end
                single_bact=single_bact2;
            end
            
        else
            single_bact=round(im); % avoid interpolated values
        end
        full_im((i-1)*grid_size+1:i*grid_size,(j-1)*grid_size+1:j*grid_size)=single_bact;
    end
end


MySample=zeros(size(full_im)+2*cell_border_offset);
MySample(cell_border_offset+1:cell_border_offset+size(full_im,1),cell_border_offset+1:cell_border_offset+size(full_im,2))=full_im;

%% Show the cell
figure(1)
clf
set(gcf,'Color','w')
imagesc(MySample);
axis image
colormap('gray')
xlabel('X [pixel]')
ylabel('Y [pixel]')
title('Virtual EColi')


% disp(['Length of bacteria [um]: ',num2str(cytoplasm_a*1e-3*pixel_size/binning)]);
% disp(['Height of bacteria [um]: ',num2str(cytoplasm_b*1e-3*pixel_size/binning)]);
% disp(['Thickness of periplasm [um]: ',num2str(periplasm_thick*1e-3*pixel_size/binning)]);

% if show==1
%     imagesc(final_im)
%     axis image
%     colormap gray
% end
% 
% if save_tiffile==1
%     imwrite(uint8(final_im),tiffile,'tif');
% end
