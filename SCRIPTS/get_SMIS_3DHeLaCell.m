function [Sample, Proj]=get_SMIS_3DHeLaCell(par)

ImageFile=fullfile(par.image_library,'Hela_Binary.tif');
if ~exist(ImageFile,"file")
    warndlg(['Image file: ',ImageFile,' not found !'])
    Sample=[];
    return
end

bg_id=0; % Id for background
cyto_id=1; % Id for cytoplasm
nuc_id=2; % Id for nucleus
membrane_id=3; % Id for plasmic membrane
nuc_membrane_id=4; % Id for plasmic membrane
cluster_id=5; % Id for membrane receptor clusters

x_dim=par.x_dim; % # of pixels in x dimension
y_dim=par.y_dim; % # of pixels in y dimension
z_dim=par.z_dim; % # of pixels in z dimension

save_projection_image=1;

plasmic_layer_thick=par.plasmic_layer_thick; % Thickness of plasmic membrane [pixel]
nuc_layer_thick=par.nuc_layer_thick; % Thickness of nuclear membrane [pixel]

n_clusters=par.n_clusters; % # of receptors clusters
cluster_diam=par.cluster_diam; % [pixels] of receptors clusters
cluster_pos_id=par.cluster_pos_id; % id of area where to position receptors clusters


cyto_layer_thick=1; % Thickness of cytoplasmic layer covering the nucleus


%%
im=imread(ImageFile);
im(im<128)=0;
im(im>=128 & im<255)=128;
im(im==128)=cyto_id;
im(im==255)=nuc_id;
im(im==0)=bg_id;

imshow(im)

Sample=zeros(x_dim,y_dim,z_dim);

%Define bottom plane
s_im=size(im);
Sample(round(x_dim-s_im(1)):round(x_dim-s_im(1))+s_im(1)-1,round(y_dim-s_im(2)):round(y_dim-s_im(2))+s_im(2)-1,1)=im;
plane=squeeze(Sample(:,:,1));
% For each cytoplasm pixel, draw a z fonction

w_bg=find(plane==0); % Indices for background
w_cyt=find(plane==cyto_id); % Indices for cytoplasm
w_nuc=find(plane==nuc_id); % Indices for nucleus
[x_bg,y_bg]=ind2sub(size(plane),w_bg);
[x_cyt,y_cyt]=ind2sub(size(plane),w_cyt);
[x_nuc,y_nuc]=ind2sub(size(plane),w_nuc);


%% Thicken the cytoplasm
disp('Creating cytoplasm ...');
% Define a steeply and then slowly rising function
% x=0:200;
x0=10;
p1=10;
p2=0.05;
% z=p1*(1-exp(-x/x0))+p2*x;
% plot(z)

MyWaitBar = waitbar(0,'Creating cytoplasm ...');

for i=1:numel(w_cyt)
    if i/100==fix(i/100)
        waitbar(i/numel(w_cyt),MyWaitBar);
    end

    [x,y]=ind2sub(size(plane),w_cyt(i));
    %Get distance to extracellular space
    d2=(x-x_bg).^2+(y-y_bg).^2;
    min_d=sqrt(min(d2));
  
    z=round(p1*(1-exp(-min_d/x0))+p2*min_d);
    Sample(x,y,1:z)=cyto_id;
end

close(MyWaitBar)

%% Thicken the nucleus
disp('Creating nucleus ...');

% Define a steeply and then slowly rising function
% x=0:200;
x0=10;
p1=20;
p2=0.1;
% z=p1*(1-exp(-x/x0))+p2*x;
% plot(z)

MyWaitBar = waitbar(0,'Creating nucleus ...');

for i=1:numel(w_nuc)
    if i/100==fix(i/100)
        waitbar(i/numel(w_nuc),MyWaitBar);
    end
    [x,y]=ind2sub(size(plane),w_nuc(i));
    %Get distance to extracellular space
    d2=(x-x_cyt).^2+(y-y_cyt).^2;
    min_d=sqrt(min(d2));
    
    w_min_d=find(d2==min(d2),1);
    z0=Sample(x_cyt(w_min_d),y_cyt(w_min_d),:);
    w_z0=find(z0>0,1,'last');
    
    z=round(p1*(1-exp(-min_d/x0))+p2*min_d);
    Sample(x,y,1:z+w_z0)=nuc_id;
    
end
close(MyWaitBar)


%% add a cytoplasmic layer above everything
disp('Adding a cytoplasmic layer ...');

w_cell=find(plane~=bg_id); % Indices for cell
[x_cell,y_cell]=ind2sub(size(plane),w_cell);

MyWaitBar = waitbar(0,'Adding a cytoplasmic layer ...');

for i=1:numel(w_cell)
    if i/100==fix(i/100)
        waitbar(i/numel(w_cell),MyWaitBar);
    end

    z0=Sample(x_cell(i),y_cell(i),:);
    w_z0=find(z0>0,1,'last');
    Sample(x_cell(i),y_cell(i),w_z0:w_z0+cyto_layer_thick)=cyto_id;
end
close(MyWaitBar)

%% add a plasmic membrane
disp('Adding a plasmic membrane...');

MyWaitBar = waitbar(0,'Adding a plasmic membrane ...');

for i=1:numel(w_cell)
    if i/100==fix(i/100)
        waitbar(i/numel(w_cell),MyWaitBar);
    end

    z0=Sample(x_cell(i),y_cell(i),:);
    w_z0=find(z0>0,1,'last');
    Sample(x_cell(i),y_cell(i),w_z0:w_z0+plasmic_layer_thick)=membrane_id; % upper membrane
    Sample(x_cell(i),y_cell(i),1:plasmic_layer_thick)=membrane_id; % lower membrane
    
  
    %Get distance to extracellular space
    d2=(x_cell(i)-x_bg).^2+(y_cell(i)-y_bg).^2;
    min_d=sqrt(min(d2));
  
    if min_d==1 % If pixel is flanking extracellular space start adding plasmic membrane on the side
        w_z0=find(z0>0,1,'last');
        w_min_d=find(d2==min(d2),1);
        Sample(x_bg(w_min_d),y_bg(w_min_d),1:w_z0)=membrane_id;
    end
    
end

close(MyWaitBar)

%% add a nuclear membrane
disp('Adding a nuclear membrane...');

MyWaitBar = waitbar(0,'Adding a nuclear membrane ...');

for i=1:numel(w_nuc)
    if i/100==fix(i/100)
        waitbar(i/numel(w_nuc),MyWaitBar);
    end

    z0=Sample(x_nuc(i),y_nuc(i),:);
    w_z0=find(z0==nuc_id,1,'last'); % Above
    Sample(x_nuc(i),y_nuc(i),w_z0-nuc_layer_thick:w_z0)=nuc_membrane_id; % Upper membrane
    w_z1=find(z0==nuc_id,1,'first'); % Below
    Sample(x_nuc(i),y_nuc(i),w_z1:w_z1+nuc_layer_thick)=nuc_membrane_id; % Lower membrane
    
  
    %Get distance to cytoplasm
    d2=(x_nuc(i)-x_cyt).^2+(y_nuc(i)-y_cyt).^2;
    min_d=sqrt(min(d2));
  
    if min_d==1 % If pixel is flanking cytoplasm space start adding nuclear membrane on the side
        z0=Sample(x_nuc(i),y_nuc(i),:);
        w_z0=find(z0==nuc_id,1,'last');
        w_z1=find(z0==nuc_id,1,'first');

        w_min_d=find(d2==min(d2),1);
        Sample(x_cyt(w_min_d),y_cyt(w_min_d),w_z1:w_z0)=nuc_membrane_id;
    end
    
end

close(MyWaitBar)


%% add small receptor clusters
% Pick random position in plasmic membrane


if n_clusters>1

    disp('Creating clusters ...');

    switch cluster_pos_id
        case 'Background'
            w_clus=find(Sample==bg_id); % Indices for background
        case 'Cytoplasm'
            w_clus=find(Sample==cyto_id); % Indices for cytoplasm
        case 'Nucleus'
            w_clus=find(Sample==nuc_id); % Indices for nucleus
        case 'Plasmic membrane'
            w_clus=find(Sample==membrane_id); % Indices for plasmic membrane
        case 'Nuclear membrane'
            w_clus=find(Sample==nuc_membrane_id); % Indices for nuclear membrane
    end
    
    [x_clus,y_clus,z_clus]=ind2sub(size(Sample),w_clus); % x, y, z for selected compartment
    r=randperm(numel(w_clus));
    w_clus_selected=w_clus(r(1:n_clusters));

    [x,y,z]=ind2sub(size(Sample),w_clus_selected);


    MyWaitBar = waitbar(0,'Adding receptor clusters ...');

    for i=1:n_clusters
        if i/100==fix(i/100)
            waitbar(i/numel(n_clusters),MyWaitBar);
        end

        d2=(x(i)-x_clus).^2+(y(i)-y_clus).^2+(z(i)-z_clus).^2;
        w=find(d2<=(cluster_diam/2)^2);

        Sample(x_clus(w),y_clus(w),z_clus(w))=cluster_id;

    end
    close(MyWaitBar)
end


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

title('Virtual 3D HeLa Cell')

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

