function MySample=get_SMIS_2D_NPCs(par)

n=par.x_dim;
m=par.y_dim;

% n=512; % image size x
% m=512; % image size y
%
N=par.N; % # NPCs randomly distributed over the field of view

l=par.nup96_size;  % [pixels] Nup96 size (area is l^2)
s_in=par.inner_radius; % [pixels] NPC inner radius
s_out=par.outer_radius; % [pixels] NPC outer radius ie diameter = 14 pixel = 110nm, ie with binning = 16, pixel size = 16/14*100 = 125 nm

% l=1;  % [pixels] Nup96 size (area is l^2)
% s_in=5; % [pixels] NPC inner radius
% s_out=7; % [pixels] NPC outer radius ie diameter = 14 pixel = 110nm, ie with binning = 16, pixel size = 16/14*100 = 125 nm

min_sep=par.npc_sep; % [pixels] Minimum separation between two NPC's
border=par.border; % border protection

% min_sep=12; % [pixels] Minimum separation between two NPC's
% border=50; % border protection

qPALM_option=par.qPALM;

na=0; % # of aggregates (simulated spurious aggregates)
a_s=20; % aggregate size

% out_image='100_NPCs_randomly_distributed.tif';
% 
% out_filedir='C:\Users\bourgeoi\Documents\MATLAB\SIMULATION\PALM\SOFTWARE\SMIS_vsn1.1\PATTERNS\2D\QPALM\';

%%

MySample=zeros(n,m);
im2=MySample; % Mask to ensure that the NPC's are well separated
sep=round(min_sep+s_out); % That's the minimum distance between two NPC centers
border=max([border,sep]); % To avoid problems at image borders

n2=n-2*border;
m2=m-2*border;

% step_x = round(n2/p);
% step_y = round(m2/q);

% X=randi(n2,N,1);
% Y=randi(m2,N,1);

%Define NUP96s positions
v=1; % pattern id, to be increased for qPALM

max_n_trials=1e+5; % In case of too high density
n_trial=1;
for i=1:N
    %     for j=1:q
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
    for theta=0:2*pi/8:(2*pi-2*pi/8)
        x=1+X + s_in*cos(theta); x=round(x);
        y=1+Y + s_in*sin(theta); y=round(y);
        MySample(x:x+l-1, y:y+l-1)=v;
        if qPALM_option==1
            v=v+1;
        end
    end
    for theta=0:2*pi/8:(2*pi-2*pi/8)
        x=1+X  + s_out*cos(theta); x=round(x);
        y=1+Y  + s_out*sin(theta); y=round(y);
        MySample(x:x+l-1, y:y+l-1)=v;
        if qPALM_option==1
            v=v+1;
        end
    end
    im2(X-sep:X+sep,Y-sep:Y+sep)=1; % Mask the corresponding region   
end

%Eventually place aggregates randomly
xy_a=rand(2,na);
x_a=n2*xy_a(1,:);
y_a=m2*xy_a(2,:);

for i=1:na
    MySample(fix(x_a -a_s):fix(x_a ),fix(y_a -a_s):fix(y_a ))=2^16-1;
end

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

