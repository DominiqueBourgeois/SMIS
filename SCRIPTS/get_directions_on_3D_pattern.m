function [u_mean, im]=get_directions_on_3D_pattern(ind, x0, y0, z0, s, show, im, D)

% PURPOSE:
% Get unit vectors from a point x0, y0, z0 to a series of indices in a 3D
% pattern, group these vectors in clusters corresponding to features in
% a pattern (ie roads in cytoskeleton structures) and get the final average
% unit vectors for these features
%
% INPUTS:
%   ind: the array of accessible indices in the image
%	x0: the x-coordinate on image [raster]
%	y0: the y-coordinates on image [raster]
%	z0: the z-coordinates on image [raster]
%	s: the image size
%	show: set to 1 to see on image
%   im: (optional] the image if show=1
%   D: Search diameter, only used if show=1
%
% OUTPUTS:
%	u_mean: average unit vectors for detected features
%   im: the updated image if show=1, otherwise []
%
% MODIFICATION HISTORY:
%	D.Bourgeois, November 2020: version > simulate_palm_vsn16.3

if nargin<7
    im=[];
    if show==1
        error('Need to input image in show=1 mode !');
    end
end

if show==1 % Set up image if we want to see it
    I=zeros(s);
    I(round(x0), round(y0), round(z0))=15;
    I(ind)=20;
end


%Get the angles into groups.
if numel(ind)>1 % Only do it if necessary
    gt=zeros(numel(ind),2); % Create the group array of angles
    gt(:,1)=ind; % Create the group array
    gt(1,2)=1; % Assign group number =1 to the first element
    
    gt2=gt; % Create a copy of gt that will be progressively shrunk over the cycles
    if show==1
        im(ind(1))=30;
    end
    
    
    s1=s(1); % The size of a line along the second dimension of the 3D image
    s2=s(1)*s(2); % The size of a plane along the third dimension of the 3D image
    
    gn=1; % Initialize Group #
    w1=0; % Initialize w1
    
    while size(gt2,1)>1 && ~isempty(w1)
        w0=find(gt2(:,2)==gn,1); % Find the first occurrence with the current group number
        while ~isempty(w0) && ~isempty(w1)
            w1=find(gt2(:,2)==0); % Restrict search to those pixels that have not been assigned a group number
            w=gt2(w0,1); % The current index
            prox_3D_coord=[w-1,w+1,w-1+s1,w+s1,w+1+s1,w-1-s1,w-s1,w+1-s1, ...
                w-1+s2,w+s2, w+1+s2,w-1+s1+s2,w+s1+s2,w+1+s1+s2,w-1-s1+s2,w-s1+s2,w+1-s1+s2, ...
                w-1-s2,w-s2, w+1-s2,w-1+s1-s2,w+s1-s2,w+1+s1-s2,w-1-s1-s2,w-s1-s2,w+1-s1-s2];
            
            w2=ismember(gt2(w1,1),prox_3D_coord);
            
            if any(w2) % Neighboring pixelx have been found
                gt2(w1(w2),2)=gn; % Assign the current group number to those pixels
                w3=ismember(gt(:,1),gt2(w1(w2),1)); % Get the corresponding indices in the original gt
                gt(w3,2)=gn; %
            end
            gt2(w0,:)=[]; % Remove the treated pixel
            w0=find(gt2(:,2)==gn,1); % Find the first occurrence with the current group number
        end
        % No neighboring pixel has been found
        if ~isempty(w1)
            gn=gn+1; % Increase group number by one
            % gt2(w1(1),2)=gn; 
            % gt(gt(:,1)==gt2(w1(1),1),2)=gn;
            w4=find(gt2(:,2)==0,1); % Pick the first pixel that has not been assign a group number
            gt2(w4,2)=gn;
            gt(gt(:,1)==gt2(w4,1),2)=gn;
        end
    end
    
    if show==1
        for i=1:gn
            im(gt(gt(:,2)==i,1))=30+10*(i-1);
            I(gt(gt(:,2)==i,1))=30+10*(i-1);
        end
    end
    
    %Get the mean directions for the different groups
    u_mean=zeros(3,gn); % The average direction and dispersion for each group
    
    for k=1:gn
        [x, y, z]=ind2sub(s,gt(gt(:,2)==k,1)); % x,y coordinates corresponding to ind
        
        %Get the Unit vectors from x,y,z
        v=[x'-round(x0); y'-round(y0); z'-round(z0)];
        u=v./vecnorm(v);
        u=u(:,~isnan(u(1,:))); %Potentially remove nan data if we are at the border of the image
        u_mean(:,k)=mean(u,2); % And get the average unit vector
    end
else
    u_mean=[];
end

if show==1 % Set up image if we want to see it
    xc=round(x0);
    yc=round(y0);
    zc=round(z0);
    %Get a local view of the image
    I_loc=I(max([1,xc-D]):min([xc+D,s(1)]),max([1,yc-D]):min([yc+D,s(2)]),max([1,zc-D]):min([zc+D,s(3)]));
    figure(1)
    clf
    vol3d('cdata', I_loc);
    alphamap('default')
    axis image
    view(3);
    set(gcf,'color','white')
    title('3D View of the selected regions')
end












