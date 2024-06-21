function crossing_found=get_potential_crossing_XYZ(xyz_s,xyz_e,sp_s,sp_e,sp_ids,w,S)

% PURPOSE:
% Check potential crossing through undesired sub patterns by a single molecule diffusing in 2D
% Useful in the case of convex shapes or for hop diffusion
% INPUTS:
%	xyz_s: the xy-coordinate on high-resolution image (Start position)
%	xyz_e: the xy-coordinate on high-resolution image (End position)
%	sp_s: sub-pattern ID at start
%	sp_e: sub-pattern ID at end
%	sp_ids: sub-pattern IDs
%   w: indices of virtual sample subpatterns
%	S: the image size
%
% OUTPUTS:
%	crossing_found: true if a crossing was found
%
% MODIFICATION HISTORY:
%	D.Bourgeois, January 2024

crossing_found=0; % Initialize

%Find the biggest slope
w_max_slope=find(abs(xyz_e-xyz_s)==max(abs(xyz_e-xyz_s)),1);

switch w_max_slope
    case 1 % max slope in X
        if xyz_e(1)>xyz_s(1)
            x1=xyz_s(1); y1=xyz_s(2); z1=xyz_s(3);
            x2=xyz_e(1); y2=xyz_e(2); z2=xyz_e(3);
        else
            x1=xyz_e(1); y1=xyz_e(2); z1=xyz_e(3);
            x2=xyz_s(1); y2=xyz_s(2); z2=xyz_s(3);
        end
        % Only proceed if there is a pixel to test in between X1 and X2
        if round(x2)>round(x1)+1
            w_ok=vertcat(w(sp_ids==sp_s | sp_ids==sp_e).w);
            % Get the affine function linking starting and ending positions
           
            % Parametric equation of the straight line xyz(i)=abc(i)*t+[x1,y1,z1];
            x=round(x1)+1:round(x2)-1; % The x values to test
            t=(x-x1)/(x2-x1);
            y=(y2-y1)*t+y1;
            z=(z2-z1)*t+z1;
           
            w_xyz=sub2ind(S,round(x),round(y),round(z));
            if any(~ismember(w_xyz,w_ok))
                crossing_found=1;
            end
        end
    case 2 % max slope in Y
        if xyz_e(2)>xyz_s(2)
            x1=xyz_s(1); y1=xyz_s(2); z1=xyz_s(3);
            x2=xyz_e(1); y2=xyz_e(2); z2=xyz_e(3);
        else
            x1=xyz_e(1); y1=xyz_e(2); z1=xyz_e(3);
            x2=xyz_s(1); y2=xyz_s(2); z2=xyz_s(3);
        end
        % Only proceed if there is a pixel to test in between X1 and X2
        if round(y2)>round(y1)+1
            w_ok=vertcat(w(sp_ids==sp_s | sp_ids==sp_e).w);
            % Get the affine function linking starting and ending positions
           
            % Parametric equation of the straight line xyz(i)=abc(i)*t+[x1,y1,z1];
            y=round(y1)+1:round(y2)-1; % The x values to test
            t=(y-y1)/(y2-y1);
            x=(x2-x1)*t+x1;
            z=(z2-z1)*t+z1;
           
            w_xyz=sub2ind(S,round(x),round(y),round(z));
            if any(~ismember(w_xyz,w_ok))
                crossing_found=1;
            end
        end
    case 3 % max slope in Z
         if xyz_e(3)>xyz_s(3)
            x1=xyz_s(1); y1=xyz_s(2); z1=xyz_s(3);
            x2=xyz_e(1); y2=xyz_e(2); z2=xyz_e(3);
        else
            x1=xyz_e(1); y1=xyz_e(2); z1=xyz_e(3);
            x2=xyz_s(1); y2=xyz_s(2); z2=xyz_s(3);
        end
        % Only proceed if there is a pixel to test in between X1 and X2
        if round(z2)>round(z1)+1
            w_ok=vertcat(w(sp_ids==sp_s | sp_ids==sp_e).w);
            % Get the affine function linking starting and ending positions
           
            % Parametric equation of the straight line xyz(i)=abc(i)*t+[x1,y1,z1];
            z=round(z1)+1:round(z2)-1; % The z values to test
            t=(z-z1)/(z2-z1);
            y=(y2-y1)*t+y1;
            x=(x2-x1)*t+x1;
           
            w_xyz=sub2ind(S,round(x),round(y),round(z));
            if any(~ismember(w_xyz,w_ok))
                crossing_found=1;
            end
        end
end
