function [x,y] = get_circle_coordinates(r)

% PURPOSE:
% Get x,y coordinates of circle of radius r
%
% INPUTS:
%	r: the circle radius [pixels]
%
% OUTPUTS:
%	x,y: the 2D cartesian coordinates [pixels]
%
% MODIFICATION HISTORY:
%	D.Bourgeois, November 2020: version > simulate_palm_vsn16.3

D=2*max([1,ceil(r)])+1; % The search diameter has to exceed what can be accessed in reality

%define x,y values around x0, y0, z0 within radius r
[x,y]=meshgrid(1:D,1:D);

% recenter
x=x-(D+1)/2;
y=y-(D+1)/2;

% keep the central circle
xy_ok=(x.^2+y.^2)<D^2/4 & (x.^2+y.^2)>(D-2)^2/4; % Empirical formula: D-2 for the lower bounderies was found by trial and error

x=x(xy_ok);
y=y(xy_ok);

