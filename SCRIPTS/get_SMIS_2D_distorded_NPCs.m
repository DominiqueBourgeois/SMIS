function MySample=get_SMIS_2D_distorded_NPCs(par)

% produces a 2D .tif file of ground truth NPCs. The NPC coordinates are
% produced by cir4mics https://github.com/uhlmanngroup/cir4mics

n=par.x_dim;
m=par.y_dim;

min_sep=par.npc_sep; % [pixels] Minimum separation between two NPC's
border=par.border; % border protection

qPALM_option=par.qPALM;

R = par.raster/par.binning; % Resolution of the produced .tif image

%% Read the NPC coordinates. All the NPCs are overlaid at the same central position
data = par.used_data;

% Extract X, Y, Z and NPC # columns
X = data.x;
Y = data.y;
% Z = data.z;
ids = data.particle+1; % Starts at 0, so add +1

N=max(ids); % # NPCs randomly distributed over the field of view

% Convert coordinates to indices based on resolution R
X_idx = round(X / R);
Y_idx = round(Y / R);
% Z_idx = round(Z / R);

%% Place the NPCs on FOV

MySample=zeros(n,m);
im2=MySample; % Mask to ensure that the NPC's are well separated
s_out=max([max(X_idx)-min(X_idx),max(Y_idx)-min(Y_idx)]); % NPC diameter in Pixels
sep=round(min_sep+s_out); % That's the minimum distance between two NPC centers
border=max([border,sep]); % To avoid problems at image borders

n2=n-2*border;
m2=m-2*border;


%Define NPC center positions
v=1; % pattern id, to be increased for qPALM
n_bad_NPCs=0; % # of bad NPCs
max_n_trials=1e+5; % In case of too high density
n_trial=1;
for i=1:N
    sep_ok=0;
    while ~sep_ok==1 && n_trial<max_n_trials
        X=randi(n2)+ border;
        Y=randi(m2)+ border;
        if im2(X,Y)==0
            sep_ok=1;
        end
        n_trial=n_trial+1;
    end

    if sep_ok~=1
        warndlg('Could not draw NPC field: reduce separation between NPCs or number of NPCs or increase image size !')
        MySample=[];
        return
    end

    x=1+X+X_idx(ids==i);
    y=1+Y+Y_idx(ids==i);
    w=sub2ind(size(MySample), x, y);
    %Check if every coordinate of the NPC could be created (if image resolution
    %is too low, several coordinates could be assigned to the same pixel)
    if numel(unique(w))~=numel(x)
        disp(['NPC #',num2str(i),' could not be created fully at this image resolution !'])
        n_bad_NPCs=n_bad_NPCs+1;
    else
        if qPALM_option==1
            MySample(w) = v:v+numel(x)-1;
            v=v+numel(x);
            % MySample(sub2ind(size(MySample), x, y)) = v; % This could be used to tag each NPC  with a single subpattern id
            % v=v+1; % This could be used to tag each NPC  with a single subpattern id
        else
            MySample(w) = 1;
        end
    end
    
    im2(X-sep:X+sep,Y-sep:Y+sep)=1; % Mask the corresponding region
end

if n_bad_NPCs>0
    MyMessage=[num2str(n_bad_NPCs),' NPCs could ne be created because of too poor image resolution !'];
    MyDlg=warndlg(MyMessage);
    waitfor(MyDlg)
end

% Make sure qPALM patterns start at 1
% if qPALM_option==1
%     u_val=unique(MySample);
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
%         newImage = zeros(size(MySample));
%         for i = 1:numel(MySample)
%             newImage(i) = valueMap(MySample(i));
%         end
%         MySample=newImage;
%     end
% end


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
title('Virtual NPCs')

