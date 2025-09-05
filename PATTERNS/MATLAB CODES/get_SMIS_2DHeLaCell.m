function MySample=get_SMIS_2DHeLaCell(par)

ImageFile=fullfile(par.image_library,'Hela_Binary.tif');
if ~exist(ImageFile,"file")
    warndlg(['Image file: ',ImageFile,' not found !'])
    MySample=[];
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
plasmic_layer_thick=par.plasmic_layer_thick; % Thickness of plasmic membrane [pixel]
nuc_layer_thick=par.nuc_layer_thick; % Thickness of nuclear membrane [pixel]

n_clusters=par.n_clusters; % # of receptors clusters
cluster_diam=par.cluster_diam; % [pixels] of receptors clusters
cluster_pos_id=par.cluster_pos_id; % id of area where to position receptors clusters

%%

image_border=round(plasmic_layer_thick/2)+3; % Define a border taking into account the thickness of the plasmic membrane
im=imresize(imread(ImageFile),[x_dim-2*image_border y_dim-2*image_border]) ;
im(im<128)=bg_id;
im(im>=128 & im<255)=cyto_id;
im(im==255)=nuc_id;

MySample=zeros(x_dim,y_dim);

%Insert the cell image into MySample
s_im=size(im);
MySample(round(image_border):round(image_border)+s_im(1)-1,round(image_border):round(image_border)+s_im(2)-1,1)=im;


%% add plasmic & nuclear membranes
disp('Adding plasmic and nuclear membranes ...');

MySample = double(MySample);

% Create binary masks for each zone
Nuc = (MySample == nuc_id);
Cyt = (MySample == cyto_id);
Bg = (MySample == bg_id);

% Define structuring elements for dilation
seNuc = strel('disk', nuc_layer_thick);
seCyt = strel('disk', plasmic_layer_thick);

% Dilate the central zone and subtract the original to get the frontier
dilatedNuc = imdilate(Nuc, seNuc);
frontierNucCyt = dilatedNuc & Cyt;

% Dilate the intermediate zone and subtract the original to get the frontier
dilatedCyt = imdilate(Cyt, seCyt);
frontierCytBg = dilatedCyt & Bg;

% Assign values to the frontiers in the original image
MySample(frontierNucCyt) = nuc_membrane_id;
MySample(frontierCytBg) = membrane_id;

% Ensure the output is an image matrix of type uint8
MySample = uint8(MySample);


%% add small receptor clusters
if n_clusters>1
    disp('Creating clusters ...');

    switch cluster_pos_id
        case 'Background'
            w_clus=find(MySample==bg_id); % Indices for background
        case 'Cytoplasm'
            w_clus=find(MySample==cyto_id); % Indices for background
        case 'Nucleus'
            w_clus=find(MySample==nuc_id); % Indices for background
        case 'Plasmic membrane'
            w_clus=find(MySample==membrane_id); % Indices for background
        case 'Nuclear membrane'
            w_clus=find(MySample==nuc_membrane_id); % Indices for background
    end

    [x_clus,y_clus]=ind2sub(size(MySample),w_clus);
    r=randperm(numel(w_clus));
    w_selected=r(1:n_clusters);

    for i=1:n_clusters
        x=x_clus(w_selected(i));
        y=y_clus(w_selected(i));
        d2=(x-x_clus).^2+(y-y_clus).^2;
        w=find(d2<=(cluster_diam/2)^2);

        MySample(x_clus(w),y_clus(w))=cluster_id;

    end
end

%%
disp('Done !');


%% Show the cell
figure(1)
clf
set(gcf,'Color','w')
imagesc(MySample);
axis image
colormap('gray')
xlabel('X [pixel]')
ylabel('Y [pixel]')
title('Virtual HeLa Cell')

