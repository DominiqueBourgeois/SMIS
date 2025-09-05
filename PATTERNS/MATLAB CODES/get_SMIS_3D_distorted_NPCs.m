function [im,Proj]=get_SMIS_3D_distorted_NPCs(par)

% produces a 3D .tif file of ground truth distorted NPCs. The NPC coordinates are
% produced by cir4mics https://github.com/uhlmanngroup/cir4mics

x_dim=par.y_dim;
y_dim=par.x_dim;
z_dim=par.z_dim;

R = par.raster/par.binning; % Resolution of the produced .tif image

qPALM_option=par.qPALM; %Set to 1 to assign a different pixel value to all nup96s in the field for qPALM experiments
min_sep=par.min_sep; % [pixels] Minimum separation between two NPC's

cell_border_offset=par.cell_border_offset; % [Pixels]

% Shift relative to optical axis
move_sample_to_coverslip=par.move_sample_to_coverslip;
z_offset=par.z_offset;

id_npc=1; % Pattern id associated to mt;

save_projection_image=1;


%% Read the NPC coordinates. All the NPCs are overlaid at the same central position
data = par.used_data;

% Extract X, Y, Z and NPC # columns
X = data.x;
Y = data.y;
Z = data.z;
ids = data.particle+1; % Starts at 0, so add +1

N=max(ids); % # NPCs randomly distributed over the field of view

% Convert coordinates to indices based on resolution R
X_idx = round(X / R);
Y_idx = round(Y / R);
Z_idx = round(Z / R);
Z_idx = Z_idx-min(Z_idx)+1;


v=1;
n_bad_NPCs=0; % # of bad NPCs


%% Full image
im3D=zeros(x_dim,y_dim,z_dim); % The full 3D image
im2D=zeros(x_dim,y_dim); % Mask to ensure that the NPC's are well separated

s_out=max([max(X_idx)-min(X_idx),max(Y_idx)-min(Y_idx)]); % NPC diameter in Pixels
sep=round(min_sep+s_out); % That's the minimum distance between two NPC centers
border=max([cell_border_offset,sep]); % To avoid problems at image borders

n2=x_dim-2*border;
m2=y_dim-2*border;

MyWaitBar = waitbar(0,'Generating NPCs ...');

% Place each individual NPC
for i=1:N
    waitbar(i/N,MyWaitBar);

    sep_ok=0;
    trial_number=1;
    while ~sep_ok==1 && trial_number<1000
        X=randi(n2)+ border;
        Y=randi(m2)+ border;
        if im2D(X,Y)==0 
            sep_ok=1;
        end
        trial_number=trial_number+1;
    end

    if sep_ok==0
        warndlg('Could not place NPCs: increase image XY size !')
        im=[];
        Proj=[];
        return
    end

    x=1+X+X_idx(ids==i);
    y=1+Y+Y_idx(ids==i);
    z=1+Z_idx(ids==i);
    w=sub2ind(size(im3D), x, y, z);
    
    %Check if every coordinate of the NPC could be created (if image resolution
    %is too low, several coordinates could be assigned to the same pixel)
    if numel(unique(w))~=numel(x)
        disp(['NPC #',num2str(i),' could not be created fully at this image resolution !'])
        n_bad_NPCs=n_bad_NPCs+1;
    else
        if qPALM_option==1
            im3D(w) = v:v+numel(x)-1;
            v=v+numel(x);
            % im3D(sub2ind(size(im3D), x, y, z)) = v; % This could be used to tag each NPC  with a single subpattern id
            % v=v+1; % This could be used to tag each NPC  with a single subpattern id
        else
            im3D(w) = 1;
        end
    end
    
    im2D(X-sep:X+sep,Y-sep:Y+sep)=1; % Mask the corresponding 2D region
end

if n_bad_NPCs>0
    MyMessage=[num2str(n_bad_NPCs),' NPCs could ne be created because of too poor image resolution !'];
    MyDlg=warndlg(MyMessage);
    waitfor(MyDlg)
end

close(MyWaitBar)

im=im3D*id_npc;

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

% Make sure qPALM patterns start at 1
% if qPALM_option==1
%     u_val=unique(Sample);
%     if (numel(u_val)-1)~=numel(X_idx)
%         MyMessage=['Number of labeled pixels created (',num2str(numel(u_val)-1), ') not equal to number of coordinates (',num2str(numel(X_idx)),') !'];
%         disp(MyMessage);
%         MyDlg=warndlg(MyMessage);
%         waitfor(MyDlg)
%         %Reorder the clusters from 1 to numel(u_val)-1
%         % Get the unique pixel values and sort them
%         disp('Reassigning pixels ...');
%         sortedValues = sort(u_val);
% 
%         % Create a mapping from the original values to the new values
%         valueMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
%         for i = 1:length(sortedValues)
%             valueMap(sortedValues(i)) = i-1;
%         end
% 
%         % Reassign the pixel values
%         newImage = zeros(size(Sample));
%         for i = 1:numel(Sample)
%             newImage(i) = valueMap(Sample(i));
%         end
%         Sample=newImage;
%     end
% end


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

title('Virtual 3D NPCs')

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


