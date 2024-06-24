function [sms, sm_par] = get_D_ini(n_dyes, sms, sm_par, im_par)

% Get initial D values for single-molecules which are not tandem
% The issue is when there are several diffusion states from a single
% sub-pattern

for i=1:n_dyes
    % First deal with non acceptor dyes
    if ~sm_par(i).is_acceptor
        K=sm_par(i).D_rate_matrix;
        confined=sm_par(i).D_confined;
        sms(i).sm=get_diff_state(sms(i).sm, K, confined);
    end
end

% If dye is an acceptor
if im_par.fluorophore_pairing_on==1
    fluorophore_pairs=im_par.fluorophore_pairs;
    n_pairs=size(fluorophore_pairs,1);
    
    
    %Go through all pairs
    for i=1:n_pairs
        dye_id_ref=fluorophore_pairs(i,1); % The reference dye
        dye_id_target=fluorophore_pairs(i,2); % The target dyes to be assigned a diffusion state
        target_id=[sms(dye_id_target).sm.id]; % the ids of the target molecules
        matched=[sms(dye_id_ref).sm.matched]; % the ids of the target molecules matched by a ref molecule
        
        % First the sm that are properly matched to a ref dye should have the same
        % diffusion state
        [~, i_matched, i2_matched]=intersect(target_id,matched); % the ids of the matched target molecules
        newVals = num2cell([sms(dye_id_ref).sm(i2_matched).diff_state]); [sms(dye_id_target).sm(i_matched).diff_state] = newVals{:};
        
        % Second the sm that are not properly matched to a ref dye should
        % be assigned a diffusion state.
        [~, i_unmatched]=setdiff(target_id,matched); % the ids of the matched target molecules
        unmatched_sm=sms(dye_id_target).sm(i_unmatched);
        
        % Treat those that had been matched to a ref dye, but the ref dye
        % is not visible (fluorescent), due to maturation
        w_unmatched_1=find([unmatched_sm.matched]==-1);
        if ~isempty(w_unmatched_1)
            K_ref=sm_par(dye_id_ref).D_rate_matrix;
            confined_ref=sm_par(dye_id_ref).D_confined;
            unmatched_sm(w_unmatched_1)=get_diff_state(unmatched_sm(w_unmatched_1), K_ref, confined_ref);
        end
        
        % Treat those that have not been matched to a ref dye, ie extra
        % target dyes
        
        w_unmatched_2=find([unmatched_sm.matched]==0);
        if ~isempty(w_unmatched_2)
            K_target=sm_par(dye_id_ref).D_rate_matrix;
            confined_target=sm_par(dye_id_ref).D_confined;
            unmatched_sm(w_unmatched_2)=get_diff_state(unmatched_sm(w_unmatched_2), K_target, confined_target);
        end
        
        newVals = num2cell([unmatched_sm.diff_state]); [sms(dye_id_target).sm(i_unmatched).diff_state] = newVals{:};
        
    end
end
    
    
    
