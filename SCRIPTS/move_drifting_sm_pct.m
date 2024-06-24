function sm=move_drifting_sm_pct(sm, im_par)

% PURPOSE:
%	Move the molecules according to set drift
%
% INPUTS:
%   sm: the single molecules
%	im_par: the imaging parameters
%
% OUTPUTS:
%   sm: the single molecules updated for drift, current position

% MODIFICATION HISTORY:
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15
%	D.Bourgeois, September 2022, use cell arrays instead of structures
%	D.Bourgeois, June 2024, use more elaborate drift from SMIS GUI. Drift
%	sequences now available from im_par

% Get the proper indices in sm
x_idx=1;
y_idx=2;
z_idx=3;

frame=im_par.current_frame;
raster=im_par.raster;

% Get current sm coordinates on detector [raster]
if im_par.simul_3D==0 % 3D mode
    [x,y,~]=get_coordinates_on_detector_pct([sm{x_idx,:}],[sm{y_idx,:}],[], im_par.binning);
else
    [x,y,z]=get_coordinates_on_detector_pct([sm{x_idx,:}],[sm{y_idx,:}],[sm{z_idx,:}], im_par.binning);
end

%Apply rotation if needed
if ~isempty(im_par.drift.dtheta)
    x0=im_par.drift.rot_x0; % Center of rotation
    y0=im_par.drift.rot_y0;
    theta=im_par.drift.dtheta(frame);

    % Apply rotation second (only in XY plane)
    r=sqrt((x-x0).^2+(y-y0).^2);
    phi=acos((x-x0)./r);
    phi2=asin((y-y0)./r);
    phi(phi2<0)=-phi(phi2<0);

    %New x,y position [raster]
    x=r.*cos(theta+phi)+x0;
    y=r.*sin(theta+phi)+y0;
end

%Then add x,y drift 
x=x + im_par.drift.dx(frame)/raster;
y=y + im_par.drift.dy(frame)/raster;

if im_par.simul_3D==1 % 3D mode
    z=z + im_par.drift.dz(frame)/raster;
    z = min(z,im_par.nz+0.49999); % Do not get out of the field of view
    z = max(z,0.50001);
end

% Update coordinates on high res image
sm(x_idx,:)=num2cell((x-0.5)*im_par.binning+0.5);
sm(y_idx,:)=num2cell((y-0.5)*im_par.binning+0.5);

if im_par.simul_3D==1 % 3D mode
    sm(z_idx,:)=num2cell((z-0.5)*im_par.binning+0.5);
end

