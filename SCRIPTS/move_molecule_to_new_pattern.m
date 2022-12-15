function P2 = move_molecule_to_new_pattern(P,w,im_par)
% PURPOSE:
%	Move a molecule to the nearest position in a new pattern
%
% INPUTS:
%	P: the coordinates of the point [pixels]
%   w: indices of the patterns where to move
%	im_par: imaging parameters
%
% OUTPUTS:
%	P2: the new coordinates of the molecule [pixels]
%
% MODIFICATION HISTORY:
%	D.Bourgeois, April 2021

%2D case
if im_par.simul_3D==0
    sz=im_par.binning*[im_par.n, im_par.m];
    [xw, yw] = ind2sub(sz,w); % the associated index
    d=(P(1)-xw).^2+(P(2)-yw).^2;
    w_min=find(d==min(d),1);
    P2=[xw(w_min),yw(w_min)];
    
elseif im_par.simul_3D==1
  sz=im_par.binning*[im_par.n, im_par.m, im_par.nz];
    [xw, yw, zw] = ind2sub(sz,w); % the associated index
    d=(P(1)-xw).^2+(P(2)-yw).^2+(P(3)-zw).^2;
    w_min=find(d==min(d),1);
    P2=[xw(w_min),yw(w_min),zw(w_min)];     
end