function [dx,dy,dz, dtheta, cdx, cdy, cdz]=calculate_drift(drift,n,simul_3D, x_pos, y_pos, raster)

% PURPOSE:
%	Calculate x, y, z drift (does not include rotation)
%
% INPUTS:
%	drift: the drift parameters
%   n: the # of frames in data collection
%   simul_3D: 1 if 3D
%   x_pos, y_pos: the coordinates where to calculate the drift 
%   raster: pixel size [nm]
%
% OUTPUTS:
%   dx,dy,dz: the drift sequence [nm]
%   dtheta: the theta drift sequence [°]
%   cdx,cdy,cdz: the cumulated drift [nm]
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2024

% Initialize cumulated drift
cdz=[];
dtheta=[];
frames=1:n; % Frame numbers

dx=drift.x1+frames/n.*drift.x2+(frames/n).^2*drift.x3+drift.x_n*randn(1,n); % X drift in [nm]
if ~isempty(drift.x_eval)
    x_fun=eval(drift.x_eval);
    x_fun=x_fun-x_fun(1); % Prevent jumps in drift at start
    if size(x_fun,1)~=size(dx,1)
        dx=dx+x_fun';
    else
        dx=dx+x_fun;
    end
end

dy=drift.y1+frames/n.*drift.y2+(frames/n).^2*drift.y3+drift.y_n*randn(1,n); % Y drift in [nm]
if ~isempty(drift.y_eval)
    y_fun=eval(drift.y_eval);
    y_fun=y_fun-y_fun(1); % Prevent jumps in drift at start
    if size(y_fun,1)~=size(dy,1)
        dy=dy+y_fun';
    else
        dy=dy+y_fun;
    end
end

if simul_3D==1 % 3D mode
    dz=drift.z1+frames/n.*drift.z2+(frames/n).^2*drift.z3+drift.z_n*randn(1,n); % Z drift in [nm]
    if ~isempty(drift.z_eval)
        z_fun=eval(drift.z_eval);
        z_fun=z_fun-z_fun(1); % Prevent jumps in drift at start
        if size(z_fun,1)~=size(dz,1)
            dz=dz+z_fun';
        else
            dz=dz+z_fun;
        end
    end
    cdz=cumsum(dz);

else
    dz=[];
end


% Evaluate rotation if any. In that case, the cumulative drifts need to be
% calculated at every frame.
if drift.rot_theta~=0 || drift.rot_n~=0 || strcmp(drift.rot_eval,'')~=1
    x0=drift.rot_x0; % [raster]
    y0=drift.rot_y0;
    dtheta=(drift.rot_theta+drift.rot_n*randn(1,n))*pi/180; % [rad]
    if ~isempty(drift.rot_eval)
        rot_fun=eval(drift.rot_eval);
        rot_fun=rot_fun-rot_fun(1); % Prevent jumps in drift at start
        dtheta=dtheta+rot_fun; % 
    end
  
    %Initialize cumulative drift
    cdx=zeros(1,n);
    cdy=zeros(1,n);

    % Apply rotation (only in XY plane) successively
    xn=x_pos; % Initial position [raster]
    yn=y_pos;    
    
    for k=1:n
        %Starting position
        xs=xn;
        ys=yn;

        %Rotation
        r=sqrt((xs-x0).^2+(ys-y0).^2);
        phi=acos((xs-x0)./r);
        phi2=asin((ys-y0)./r);
        phi(phi2<0)=-phi(phi2<0);

        %New position
        xn=r.*cos(dtheta(k)+phi)+x0+dx(k)/raster; % [raster] (dx in [nm])
        yn=r.*sin(dtheta(k)+phi)+y0+dy(k)/raster;

        %Cumulant drift [nm]
        cdx(k)=(xn-x_pos)*raster;
        cdy(k)=(yn-y_pos)*raster;
    end

else % Just apply cumulative sum
    cdx=cumsum(dx);
    cdy=cumsum(dy);
end

