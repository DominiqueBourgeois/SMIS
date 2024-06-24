% Script to transform image data into .csv file, by clicking on image at
% defined positions.
% Useful to transform image of spectrum into real spectrum

im_file='C:\Users\bourgeoi\Documents\MATLAB\SIMULATION\PALM\SOFTWARE\SMIS_vsn2.1\FLUOROPHORES\SPECTRA\CF660C\Emission.png'; % Name of image
out_csvfile='C:\Users\bourgeoi\Documents\MATLAB\SIMULATION\PALM\SOFTWARE\SMIS_vsn2.1\FLUOROPHORES\SPECTRA\CF660C\Emission.csv'; % Name of image

x_lim=[635,805]; % X range covered: data points should be recorded at these limits
y_lim=[0,100]; % Y range covered in spectrum: data points should be recorded at these limits
reverse_y=1; % Set to 1 to reverse Y axis
reverse_x=0; % Set to 1 to reverse X axis

smooth_spectrum=1; % Set to 1 to smooth spectrum with spline fit
x_fit=x_lim(1):1:x_lim(2); % Target X coordinates

% Load an image into Matlab
close all

f=figure;
pic = imread(im_file);
imshow(pic);
clearvars a
i=1;
while true
    shg
    dcm_obj = datacursormode(1);
    set(dcm_obj,'DisplayStyle','window',...
        'SnapToDataVertex','off','Enable','on')
    waitforbuttonpress
    tmp_c=getCursorInfo(dcm_obj)
    data(i)=tmp_c;
    i=i+1;
    
    click_type=get(f,'SelectionType');
    if strcmp(click_type,'alt') %right click
        break
    end  
end

% Transform into matrix
xy = cell2mat({data.Position}');
% Sort and remove duplicates
[~,idu] = unique(xy(:,1));
xy_2 = xy(idu,:);

if reverse_y==1
    xy_2(:,2)=max(xy_2(:,2))-xy_2(:,2);
end
if reverse_x==1
    xy_2(:,1)=max(xy_2(:,1))-xy_2(:,1);
end

X_min=min(xy_2(:,1));
X_max=max(xy_2(:,1));
Y_min=min(xy_2(:,2));
Y_max=max(xy_2(:,2));


xy_2(:,1)=x_lim(1)+(xy_2(:,1)-X_min)*(x_lim(2)-x_lim(1))/(X_max-X_min);
xy_2(:,2)=y_lim(1)+(xy_2(:,2)-Y_min)*(y_lim(2)-y_lim(1))/(Y_max-Y_min);

x_lim=[280,800]; % X range covered: data points should be recorded at these limits
y_lim=[0,5e+4]; % Y range covered in spectrum: data points should be recorded at these limits

%Do spline fitting
% p = 0.001;
p=0.01;
% p = 0.5;
y_fit = csaps(xy_2(:,1),xy_2(:,2),p,x_fit);

%Interpolate
% x_interp=650:800;
% y_interp=interp1(xy_2(:,1),xy_2(:,2),x_interp,'pchip');
% y_interp=interp1(xy_2(:,1),xy_2(:,2),x_interp,'spline');

figure(1)
clf
plot(xy_2(:,1),xy_2(:,2),'red')
hold on
plot(x_fit,y_fit,'green')
% plot(x_interp,y_interp,'blue')


out_data=[x_fit; y_fit];

%Write the data
csvwrite(out_csvfile,out_data');


