function sms=move_drifting_sm(n_dyes,sms, im_par)

% PURPOSE:
%	Move the molecules according to set drift
%
% INPUTS:
%   n_dyes: the # of dyes
%   sms: the single molecules (with coordinates on high-resolution image) in raster units
%	im_par: the imaging parameters
%
% OUTPUTS:
%   sms: the single molecules updated for drift, current position

% MODIFICATION HISTORY:
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15

%
% Apply drift first
drift=im_par.xyz_drift;
frame=im_par.current_frame;
raster=im_par.raster;

dx=drift(1,1)+frame/im_par.n_images*drift(1,2)+(frame/im_par.n_images)^2*drift(1,3);
dx=(dx*(1+drift(1,4)*randn(1)))/raster; % X drift in raster unit

dy=drift(2,1)+frame/im_par.n_images*drift(2,2)+(frame/im_par.n_images)^2*drift(2,3);
dy=(dy*(1+drift(2,4)*randn(1)))/raster; % Y drift in nm

x0=im_par.rot_drift(1);
y0=im_par.rot_drift(2);
theta=im_par.rot_drift(3)*pi/180;
theta=theta*(1+im_par.rot_drift(4)*randn(1));

if im_par.simul_3D==1 % 3D mode
    %Only apply drift
    dz=drift(3,1)+frame/im_par.n_images*drift(3,2)+(frame/im_par.n_images)^2*drift(3,3);
    dz=(dz*(1+drift(3,4)*randn(1)))/raster; % Z drift in nm
end


for i=1:n_dyes
    if im_par.simul_3D==0 % 3D mode
        [x,y,~]=get_coordinates_on_detector(sms(i).sm, im_par);
    else
        [x,y,z]=get_coordinates_on_detector(sms(i).sm, im_par);
    end
    x2=x + dx;
    y2=y + dy;
    
    if im_par.simul_3D==1 % 3D mode
        z2=z + dz;
        z2 = min(z2,im_par.nz+0.49999); % Do not get out of the field of view
        z2 = max(z2,0.50001);
    end
    
    % Apply rotation second (only in XY plane)
    r=sqrt((x2-x0).^2+(y2-y0).^2);
    phi=acos((x2-x0)./r);
    phi2=asin((y2-y0)./r);
    phi(phi2<0)=-phi(phi2<0);
    
    x2=r.*cos(theta+phi)+x0;
    y2=r.*sin(theta+phi)+y0;
    
    % Update coordinates on high res image
    newVals = num2cell((x2-0.5)*im_par.binning+0.5); [sms(i).sm.x] = newVals{:};
    newVals = num2cell((y2-0.5)*im_par.binning+0.5); [sms(i).sm.y] = newVals{:};
    
    if im_par.simul_3D==1 % 3D mode
        newVals = num2cell((z2-0.5)*im_par.binning+0.5); [sms(i).sm.z] = newVals{:};    
    end
end

end
