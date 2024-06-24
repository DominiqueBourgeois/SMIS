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
im=imread(ImageFile);
im(im<128)=0;
im(im>=128 & im<255)=128;
im(im==128)=cyto_id;
im(im==255)=nuc_id;
im(im==0)=bg_id;

MySample=zeros(x_dim,y_dim);

%Define bottom plane
s_im=size(im);
MySample(round(x_dim-s_im(1)):round(x_dim-s_im(1))+s_im(1)-1,round(y_dim-s_im(2)):round(y_dim-s_im(2))+s_im(2)-1,1)=im;


w_bg=find(MySample==0); % Indices for background
w_cyt=find(MySample==cyto_id); % Indices for cytoplasm
w_nuc=find(MySample==nuc_id); % Indices for nucleus
[x_bg,y_bg]=ind2sub(size(MySample),w_bg);
[x_nuc,y_nuc]=ind2sub(size(MySample),w_nuc);

%% add a plasmic membrane
disp('Adding a plasmic membrane...');
r=(plasmic_layer_thick-1)/2;

MyWaitBar = waitbar(0,'Adding a plasmic membrane ...');

for i=1:numel(w_cyt)
    if i/100==fix(i/100)
        waitbar(i/numel(w_cyt),MyWaitBar);
    end
    [x,y]=ind2sub(size(MySample),w_cyt(i));
    %Get distance to extracellular space
    d2=(x-x_bg).^2+(y-y_bg).^2;
    min_d=sqrt(min(d2));

    w_min_d=find(d2==min(d2),1);

    if min_d==1 % We are at the border
        % Add membrane
        MySample(round(x_bg(w_min_d)-r):round(x_bg(w_min_d)+r),...
            round(y_bg(w_min_d)-r):round(y_bg(w_min_d)+r))=membrane_id;
    end
end
close(MyWaitBar)

%% add a nuclear membrane
disp('Adding a nuclear membrane...');
r=(nuc_layer_thick-1)/2;

MyWaitBar = waitbar(0,'Adding a nuclear membrane ...');

for i=1:numel(w_cyt)
    if i/100==fix(i/100)
        waitbar(i/numel(w_cyt),MyWaitBar);
    end
    [x,y]=ind2sub(size(MySample),w_cyt(i));
    %Get distance to extracellular space
    d2=(x-x_nuc).^2+(y-y_nuc).^2;
    min_d=sqrt(min(d2));

    w_min_d=find(d2==min(d2),1);

    if min_d==1 % We are at the border
        % Add nuc membrane
        MySample(round(x_nuc(w_min_d)-r):round(x_nuc(w_min_d)+r),...
            round(y_nuc(w_min_d)-r):round(y_nuc(w_min_d)+r))=nuc_membrane_id;
    end

end

close(MyWaitBar)

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

    MyWaitBar = waitbar(0,'Adding receptor clusters ...');

    for i=1:n_clusters
        if i/100==fix(i/100)
            waitbar(i/numel(n_clusters),MyWaitBar);
        end

        x=x_clus(w_selected(i));
        y=y_clus(w_selected(i));
        d2=(x-x_clus).^2+(y-y_clus).^2;
        w=find(d2<=(cluster_diam/2)^2);

        MySample(x_clus(w),y_clus(w))=cluster_id;

    end
    close(MyWaitBar)

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

