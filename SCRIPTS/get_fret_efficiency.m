function f = get_fret_efficiency(xyz1, xyz2, sm_par, par)
% NAME:
%   get_fret_efficiency
%
% PURPOSE:
%	Get fret efficiency between a donor and an acceptor molecule, knowing
%	R0
%
% CATEGORY:
%	Single molecules simulation.
%
% CALLING SEQUENCE:
%   f = get_fret_efficiency(xy1, xy2, sm_par, par);
%
% INPUTS:
%   xyz1: the [x,y,z] coordinates of the donor
%   xyz2: the [x,y,z] coordinates of the acceptor
%   sm_par: single molecule parameters containing the fieldname 'R0' 
%   par: general parameters containing the fieldname 'raster' 
%
% OUTPUTS:
%	f: the fret efficiency
%
% COMMON BLOCKS:
%	None.
%
% SIDE EFFECTS:
%	None.
%
% RESTRICTIONS:
%	None.
%
% MODIFICATION HISTORY:
%	D.Bourgeois, April 2011.


x1=xyz1(1);
y1=xyz1(2);
z1=xyz1(3);
x2=xyz2(1);
y2=xyz2(2);
z2=xyz2(3);


% compute the distance between the two molecules in A
R = 10*par.raster*sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2); % par.raster in nm

f = 1/(1+(R/sm_par.R0_D)^6);
 
if par.debug
    disp(['Distance between the 2 molecules is [A]: ',num2str(R)]);
    disp(['Fret efficiency: ',num2str(f)]);
end
end
