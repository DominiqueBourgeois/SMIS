function stop_data_collection=check_stop_pct(sms, sm_par, im_par)
%
% PURPOSE:
%   Check if data collection to be continued. Stop if all molecules
%   bleached or gone out of the field of view.
%
% INPUTS:
%   sms: The single molecules
%   sm_par: the single molecule parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.
%	D.Bourgeois, September 2022, optimized for parallel computing

%Stop data collection if all dyes have been bleached.
n_fluorophores=size(sm_par,2);
stop_id=0;
for k=1:n_fluorophores
    if sm_par(k).photophysical_status.cum_n_bleached(im_par.current_frame)==sm_par(k).n_mol_eff
        stop_id=1+stop_id;
    end
end
if stop_id==n_fluorophores
    disp('All molecules have been bleached, stopping data collection !');
    stop_data_collection=1;
    return
end

% Stop data collection if all dyes have gone out of the field of view
if im_par.add_drift==1 || im_par.add_diffusion==1
    stop_id2=0;
    x_idx=1;
    y_idx=2;
    z_idx=3;
    for k=1:n_fluorophores
        [x,y,~]=get_coordinates_on_detector_pct([sms(k).sm_cell{x_idx,:}],[sms(k).sm_cell{y_idx,:}],[sms(k).sm_cell{z_idx,:}], im_par.binning);
        n_in_fov=check_in_fov(x,y,1.2,im_par);
        if n_in_fov==0
            stop_id2=1+stop_id2;
        end
    end
    % Stop data collection if all dyes have been bleached
    if stop_id2==n_fluorophores
        disp('All molecules disappeared from the Field of View due to drift or diffusion !');
        disp('Stoping data collection !');
        stop_data_collection=1;
        return
    end
    
end

stop_data_collection=0; % If everything is okay and data collection can be continued


