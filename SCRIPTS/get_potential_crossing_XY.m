function crossing_found=get_potential_crossing_XY(xy_s,xy_e,sp_s,sp_e,sp_ids,w,S)

% PURPOSE:
% Check potential crossing through undesired sub patterns by a single molecule diffusing in 2D
% Useful in the case of convex shapes or for hop diffusion
% INPUTS:
%	xy_s: the xy-coordinate on high-resolution image (Start position)
%	xy_e: the xy-coordinate on high-resolution image (End position)
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

if abs(xy_e(1)-xy_s(1))>=abs(xy_e(2)-xy_s(2)) % The segment is more horizontal, scan horizontally
    if xy_e(1)>xy_s(1)
        x1=xy_s(1); y1=xy_s(2);
        x2=xy_e(1); y2=xy_e(2);
    else
        x1=xy_e(1); y1=xy_e(2);
        x2=xy_s(1); y2=xy_s(2);
    end
    % Only proceed if there is a pixel to test in between X1 and X2
    if round(x2)>round(x1)+1
        w_ok=vertcat(w(sp_ids==sp_s | sp_ids==sp_e).w);
        % Get the affine function linking starting and ending positions
        a=(y2-y1)/(x2-x1);
        b=y1-a*x1;

        % for x=round(x1)+1:round(x2)-1 % This is slower
        %     y=a*x+b;
        %     w_xy=sub2ind(S,round(x),round(y));
        %     if ~ismember(w_xy,w_ok)
        %         crossing_found=1;
        %         w_crossing=w_xy;
        %         return
        %     end
        % end
        x=round(x1)+1:round(x2)-1; % This is faster
        y=a*x+b;
        w_xy=sub2ind(S,round(x),round(y));
        if any(~ismember(w_xy,w_ok))
            crossing_found=1;
        end

    end

else % Scan vertically
    if xy_e(2)>xy_s(2)
        x1=xy_s(1); y1=xy_s(2);
        x2=xy_e(1); y2=xy_e(2);
    else
        x1=xy_e(1); y1=xy_e(2);
        x2=xy_s(1); y2=xy_s(2);
    end
    % Only proceed if there is a pixel to test in between Y1 and Y2
    if round(y2)>round(y1)+1
        w_ok=vertcat(w(sp_ids==sp_s | sp_ids==sp_e).w);
        % Get the affine function linking starting and ending positions
        a=(x2-x1)/(y2-y1);
        b=x1-a*y1;
        y=round(y1)+1:round(y2)-1;
        x=a*y+b;
        w_xy=sub2ind(S,round(x),round(y));
        if any(~ismember(w_xy,w_ok))
            crossing_found=1;
        end
    end
end
