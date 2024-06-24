function  [P_out, V_out]=get_p_v_3D_out(P_in, V_in, S, v, dt, raster, binning, margin_factor, w, sz)

% PURPOSE:
%   Select velocity vector which is closest to the initial direction, and that allows for the molecule to
%   stay within the pattern at its current speed. 3D case
%
% INPUTS:
%	P_in: the initial position vector [raster]
%	V_in: the initial velocity vector [raster.s-1]
%	S: the sphere to explore around current position
%	v: Speed up the molecule
%   dt: Moving time
%	raster: = im_par.raster
%	binning: = im_par.binning
%	margin_factor: = sm_par.margin_factor
%	w: = sm_pattern_indices.w_patterns(c_ds).w
%   sz: image size
%
% OUTPUTS:
%	P_out: the new position vector [raster]
%	V_out: the new velocity vector [raster.s-1]
%
% MODIFICATION HISTORY:
%	D.Bourgeois, November 2020: version > simulate_palm_vsn16.3

% Get the coordinates of the V_circle around current point
x2=P_in(1)+S.x;
y2=P_in(2)+S.y;
z2=P_in(3)+S.z;

% make sure circle is inside image
x2(x2<1)=1+margin_factor;
x2(x2>sz(1))=sz(1)-margin_factor;
y2(y2<1)=1+margin_factor;
y2(y2>sz(2))=sz(2)-margin_factor;
z2(z2<1)=1+margin_factor;
z2(z2>sz(3))=sz(3)-margin_factor;

%look at the intersection
int=intersect(sub2ind(sz,round(x2),round(y2),round(z2)),w); % the indices of the circle on the current pattern

if ~isempty(int) % This should be the case unless the pattern is too small locally
    %The possible reached pixels are
    [tp_x, tp_y, tp_z]=ind2sub(sz,int); % x,y, z coordinates corresponding to ind_cp_s
    
    %The possible new vectors on the unit circle are
    V=[tp_x'-round(P_in(1)); tp_y'-round(P_in(2)); tp_z'-round(P_in(3))];
    U=V./vecnorm(V);
    
    %The angles between current and new possible velocity vectors are
    dth=acos(dot(repmat(V_in'/norm(V_in),[1,size(U,2)]),U));
    
    %We select the smallest abs(dth) available
    w_th=find(abs(dth)==min(abs(dth)),1); % The chosen theta is the first one for which the cumulative sum exceeds a random number between 0 an 1
    
    %Redefine the velocity vector
    V_out=U(:,w_th)'*1e+3*v/raster*binning; % in [raster.s-1]
    
    %And update the position
    P_out=P_in+dt*V_out;
    
    %But maybe again the molecule will have gone out of the field of view
    if any([P_out(1)<1 P_out(2)<1 P_out(3)<1]) || any([P_out(1)>sz(1) P_out(2)>sz(2) P_out(3)>sz(3)]) % check if molecule is inside the FOV
        disp('Molecule may have hit image border ! ');
        warning('Setting speed to opposite current speed for current molecule !');
        V_out=-V_in; % Simply reverse the speed
        P_out=P_in; % Simply do not move the molecule
    end
else % Treat the case where the pattern may be locally two small, then the molecule should be assigned a velocity of zero
    disp(['current pattern does not allow molecule to move at speed: ', num2str(v), ' um/s'])
    warning('Setting speed to zero !');
    V_out=[0,0,0]; % Simply reverse the speed
    P_out=P_in; % Simply do not move the molecule
end
