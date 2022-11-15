function MySample=get_SMIS_2DLinePattern(par)


n=par.x_dim; %number of column
m=par.y_dim; %number of slices

nz=par.z_dim; % Only used for 3D simulations; must be >3 if create_pattern=1;

%the image will consist of a set of lines
%spacing between lines [pixel]
h_line_space=par.h_line_space;
v_line_space=par.v_line_space;

%thickness of lines [pixel]
% line_thick=1000;
h_line_thick=par.h_line_thick;
v_line_thick=par.v_line_thick;

%position of first line from top [nm]
% first_line_pos=1;
% first_line_pos=par.first_line_pos;

%number of lines to draw
% number_of_lines=26;
h_number_of_lines=par.h_number_of_lines;
v_number_of_lines=par.v_number_of_lines;


%Tilted lines
include_tilted_lines=par.include_tilted_lines;
N=par.t_number_of_lines;
t_line_thick=par.t_line_thick;

%depth of pattern
% pattern_depth=[-200, 0, 200]; % [pixel] only used if simul_3D=1
pattern_depth=par.pattern_depth;

%safety border
border=par.border;

simul_3D=par.simul_3D;

line_id=1;
%%

if simul_3D==0
    MySample = double(zeros(m,n));
elseif simul_3D==1
    if fix(abs(pattern_depth)/im_par.raster)>=nz/2
        error('Increase nz to reach this pattern depth or decrease pattern_depth !');
    end
    MySample = double(zeros(m,n,nz));
end

%Figure out the pattern
% Horizontal lines
if h_number_of_lines>0
    h_thick=h_number_of_lines*h_line_thick+(h_number_of_lines-1)*h_line_space;
    h_first_line_pos=max([1,round((m-h_thick)/2)]);
    h_last_line_pos=min([m,round((m+h_thick)/2)]);

    h_line_increment=h_line_thick+h_line_space;

    for i=1:h_line_thick
        if simul_3D==0
            MySample(h_first_line_pos+i-1:h_line_increment:h_last_line_pos+h_line_thick-1, : )=line_id;
        elseif simul_3D==1
            MySample(h_first_line_pos+i-1:h_line_increment:h_last_line_pos+h_line_thick-1, : ,nz/2+pattern_depth)=line_id;
        end
    end
end

% Vertical lines
if v_number_of_lines>0
    v_thick=v_number_of_lines*v_line_thick+(v_number_of_lines-1)*v_line_space;
    v_first_line_pos=max([1,round((n-v_thick)/2)]);
    v_last_line_pos=min([n,round((n+v_thick)/2)]);

    v_line_increment=v_line_thick+v_line_space;

    for i=1:v_line_thick
        if simul_3D==0
            MySample(:,v_first_line_pos+i-1:v_line_increment:v_last_line_pos+v_line_thick-1)=1;
        elseif simul_3D==1
            MySample(:,v_first_line_pos+i-1:v_line_increment:v_last_line_pos+v_line_thick-1, nz/2+pattern_depth)=1;
        end
    end
end

%Tilted lines
if include_tilted_lines==1
    MySample_t = double(zeros(m,n));
    rot_val_XY=180; % [Deg]

    for k=1:N
        MySample_t(:)=0;
        disp(['Computing Line #: ', num2str(k)]);
        if rand<0.5
            x=rand*m;
            MySample_t(max(round(x),1):min(round(x+t_line_thick-1),m),:)=100;
        else
            y=rand*n;
            MySample_t(:,max(round(y),1):min(round(y+t_line_thick-1),n),:)=100;
        end

        r=rot_val_XY*rand; % angle between cylinder axis and x,y planeclf
        MySample=MySample | imrotate(MySample_t,r,'bilinear','crop'); % rotation around z
    end
end

MySample=MySample*line_id;

if border>0
    MySample(1:border,:)=0;
    MySample(end-border:end,:)=0;
    MySample(:,1:border)=0;
    MySample(:,end-border:end)=0;
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
title('Crossing Lines')

