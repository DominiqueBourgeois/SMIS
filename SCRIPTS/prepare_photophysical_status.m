function sm_par=prepare_photophysical_status(sm_par,n_images)


% Added sm_par.initial_fluo_state July 2020
% Modified sm_par.initial_fluo_state July 2021: was not consistent and not
% working for PAFPs

n_dyes=size(sm_par,2);

photophysical_status=struct(...
    'n_bleached',zeros(1,n_images),...
    'n_activated',zeros(1,n_images),...
    'n_blinked',zeros(1,n_images),...
    'cum_n_bleached',zeros(1,n_images),...
    'cum_n_activated',zeros(1,n_images),...
    'cum_n_blinked',zeros(1,n_images)...
    );

for k=1:n_dyes
    sm_par(k).photophysical_status=photophysical_status;
    if sm_par(k).initial_state<sm_par(k).converted_state
        sm_par(k).initial_fluo_state=min(sm_par(k).fluorescent_states);
    else
        sm_par(k).initial_fluo_state=min(sm_par(k).fluorescent_states(sm_par(k).fluorescent_states>=sm_par(k).converted_state));
    end
    
%     sm_par(k).initial_fluo_state=max(sm_par(k).fluorescent_states(sm_par(k).fluorescent_states<=sm_par(k).initial_state));
end