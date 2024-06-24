function f = search_close_subpatterns(P,Vdt,w,r,im_par)
% PURPOSE:
%	Look if a 2D or 3D subpattern defined by indices w is within r pixels from a
%	point of coordinate x,y,z
%
% INPUTS:
%	P: the coordinates of the point [pixels]
%	Vdt: the vector distance possibly reached due to velocity of the point [pixels]
%   w: indices of the patterns
%   r: search radius [pixels]
%	im_par: imaging parameters
%
% OUTPUTS:
%   f = 0 (pattern is not within r) or 1 (pattern is within r)
%
% MODIFICATION HISTORY:
%	D.Bourgeois, December 2019: version > 15.3

%2D case
if im_par.simul_3D==0
    % R=max([1,round(r)]): The search radius [pix]; has to be at least 1 to allow searching around current pixel
    D=2*max([1,round(r)])+1; % The search diameter
    
    %define x,y values around x0, y0, z0 within radius r
    [x,y]=meshgrid(1:D,1:D);
    
    % recenter
    x=x-(D+1)/2;
    y=y-(D+1)/2;
    
    % remove corners and keep the central disk
    xy_ok=(x.^2+y.^2)<D^2/4;
    
    x=x(xy_ok);
    y=y(xy_ok);
    
    % place on main image
    if Vdt<1 % No velocity or no move by more than 1 pixel
        x=P(1)+x;
        y=P(2)+y;
    else
        % In that case get the intersection of the current velocity vector with the new pattern
        n_pos=round(norm(Vdt));
        P=P+(0:n_pos)'.*Vdt/n_pos; % Possible positions on the new pattern
        x=repmat(P(:,1)',[size(x),1])+repmat(x,[1,size(P,1)]);
        y=repmat(P(:,2)',[size(y),1])+repmat(y,[1,size(P,1)]);      
    end
    
    % make sure disk is inside image
    sz=im_par.binning*[im_par.n, im_par.m];
    xy_ok=x>=1 & x<=sz(1) & y>=1 & y<=sz(2);
    x=x(xy_ok);
    y=y(xy_ok);
    
    %look at the intersection
    ind = sub2ind(sz,round(x),round(y)); % the associated index
    
    if isempty(ind)
        f=0;
    else
        f=~isempty(intersect(ind,w(w>=min(ind) & w<=max(ind))));
    end
    
elseif im_par.simul_3D==1
    D_XY=2*max([1,round(r(1))])+1; % The search diameter in XY
%     D_Z=2*max([1,round(r(2))])+1; % The search diameter in Z
    D_Z=D_XY; % The search diameter in Z
    
    if P(3)==0; error('Z-coordinate must be >= 0 !'); end
    
    [x,y,z]=meshgrid(1:D_XY,1:D_XY,1:D_Z);
    
    % recenter
    x=x-(D_XY+1)/2;
    y=y-(D_XY+1)/2;
    z=z-(D_Z+1)/2;
    
    % remove corners and keep the central ellipsoid
    
    xyz_ok=(x.^2/(D_XY/2)^2+y.^2/(D_XY/2)^2+z.^2/(D_Z/2)^2)<1;
    
    x=x(xyz_ok);
    y=y(xyz_ok);
    z=z(xyz_ok);
    
    % place on main image
    if Vdt<1 % No velocity or no move by more than 1 pixel
        x=P(1)+x;
        y=P(2)+y;
        z=P(3)+z;
    else
        % In that case get the intersection of the current velocity vector with the new pattern
        n_pos=round(norm(Vdt));
        P=P+(0:n_pos)'.*Vdt/n_pos; % Possible positions on the new pattern
        x=repmat(P(:,1)',[size(x),1])+repmat(x,[1,size(P,1)]);
        y=repmat(P(:,2)',[size(y),1])+repmat(y,[1,size(P,1)]);
        z=repmat(P(:,3)',[size(z),1])+repmat(z,[1,size(P,1)]);
    end
    
    % make sure pixels are inside image
    sz=im_par.binning*[im_par.n, im_par.m, im_par.nz];
    xyz_ok= x>=1 & x<=sz(1) & y>=1 & y<=sz(2) & z>=1 & z<=sz(3);
    x=x(xyz_ok);
    y=y(xyz_ok);
    z=z(xyz_ok);
    
    xyz=unique([round(x), round(y), round(z)],'rows');
    
    %look at the intersection
    ind = sub2ind(sz,xyz(:,1),xyz(:,2),xyz(:,3)); % the associated index
    
    if isempty(ind)
        f=0;
    else
        f=~isempty(intersect(ind,w));
    end
    
else
end