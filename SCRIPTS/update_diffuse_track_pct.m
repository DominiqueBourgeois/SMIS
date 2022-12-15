function sms=update_diffuse_track_pct(n_fluorphores, sms, sm_par, im_par)

%
% PURPOSE:
%   % Correct for xyz_track positions to the mean position instead of the end position during frame time
%
% INPUTS:
%   n_fluorphores: the # of fluorophores
%   sms: the single molecules
%	sm_par: the sm parameters
%   im_par: the imaging parameters
%
% OUTPUTS:
%   sms: the updated single molecules
%
% MODIFICATION HISTORY:
%	D.Bourgeois, January 2022.
 
for k=1:n_fluorphores

    if im_par.simul_3D==0
        w_idx = find(matches(sm_par(k).sm_fn,{'x_track','y_track','sub_x','sub_y'})==1);
        x_track_idx=w_idx(1);
        y_track_idx=w_idx(2);
        sub_x_idx=w_idx(3);
        sub_y_idx=w_idx(4);
    else
        w_idx = find(matches(sm_par(k).sm_fn,{'x_track','y_track','z_track','sub_x','sub_y','sub_z'})==1);
        x_track_idx=w_idx(1);
        y_track_idx=w_idx(2);
        z_track_idx=w_idx(3);
        sub_x_idx=w_idx(4);
        sub_y_idx=w_idx(5);
        sub_z_idx=w_idx(6);
    end

    for i=1:size(sms(k).sm_cell,2)
        if ~isempty(sms(k).sm_cell{sub_x_idx,i})
            
            [sms(k).sm_cell{x_track_idx,i}(end,1)]=mean(sms(k).sm_cell{sub_x_idx,i});
            [sms(k).sm_cell{y_track_idx,i}(end,1)]=mean(sms(k).sm_cell{sub_y_idx,i});
            
            if im_par.simul_3D==1
                [sms(k).sm_cell{z_track_idx,i}(end,1)]=mean(sms(k).sm_cell{sub_z_idx,i});
            end
        end
    end
end


