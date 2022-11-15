function sms=update_diffuse_track(sms,im_par)

%
% PURPOSE:
%   % Correct for xyz_track positions to the mean position instead of the end position during frame time
%
% INPUTS:
%   sms: the single molecules
%   im_par: the imaging parameters
%
% OUTPUTS:
%   sms: the updated single molecules
%
% MODIFICATION HISTORY:
%	D.Bourgeois, January 2022.

n_fluorphores=numel(sms);
for k=1:n_fluorphores
    for i=1:numel(sms(k).sm)
        
        if ~isempty(sms(k).sm(i).sub_x)
            
            [sms(k).sm(i).x_track(end,1)]=mean(sms(k).sm(i).sub_x);
            [sms(k).sm(i).y_track(end,1)]=mean(sms(k).sm(i).sub_y);
            
            if im_par.simul_3D==1
                [sms(k).sm(i).z_track(end,1)]=mean(sms(k).sm(i).sub_z);
            end
        end
    end
end


