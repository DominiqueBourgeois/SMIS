function [sms, sm_par] = get_V_ini(n_dyes, sms, sm_par, im_par)

% Get initial velocity vectors for single-molecules which are not tandem

if im_par.simul_3D==0
    for i=1:n_dyes
        % First deal with non acceptor dyes
        if ~sm_par(i).is_acceptor
            % eventually initialize velocities
            if any(sm_par(i).V)
                sm_par(i) = get_velocity_radii(im_par, sm_par(i));
                x_h=[sms(i).sm.x];
                y_h=[sms(i).sm.y];
                ds=[sms(i).sm.diff_state];
                sp=[sms(i).sm.c_sp];
                v_h = get_initial_directions(x_h, y_h, 0, sp, ds, im_par, sm_par(i));
                newVals=num2cell(v_h(:,1)); [sms(i).sm.v_x]=newVals{:};
                newVals=num2cell(v_h(:,2)); [sms(i).sm.v_y]=newVals{:};
            end
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
            % velocity
            [~, i_matched, i2_matched]=intersect(target_id,matched); % the ids of the matched target molecules
            newVals = num2cell([sms(dye_id_ref).sm(i2_matched).v_x]); [sms(dye_id_target).sm(i_matched).v_x] = newVals{:};
            newVals = num2cell([sms(dye_id_ref).sm(i2_matched).v_y]); [sms(dye_id_target).sm(i_matched).v_y] = newVals{:};
            
            % Second the sm that are not properly matched to a ref dye should
            % be assigned a diffusion state.
            [~, i_unmatched]=setdiff(target_id,matched); % the ids of the matched target molecules
            unmatched_sm=sms(dye_id_target).sm(i_unmatched);
            
            % Treat those that had been matched to a ref dye, but the ref dye
            % is not visible (fluorescent), due to maturation
            w_unmatched_1=find([unmatched_sm.matched]==-1);
            if ~isempty(w_unmatched_1)
                if any(sm_par(dye_id_ref).V)
                    sm_par(dye_id_ref) = get_velocity_radii(im_par, sm_par(dye_id_ref));
                    x_h=[unmatched_sm(w_unmatched_1).x];
                    y_h=[unmatched_sm(w_unmatched_1).y];
                    ds=[unmatched_sm(w_unmatched_1).diff_state];
                    sp=[unmatched_sm(w_unmatched_1).c_sp];
                    v_h = get_initial_directions(x_h, y_h, 0, sp, ds, im_par, sm_par(dye_id_ref));
                    newVals=num2cell(v_h(:,1)); [unmatched_sm(w_unmatched_1).v_x]=newVals{:};
                    newVals=num2cell(v_h(:,2)); [unmatched_sm(w_unmatched_1).v_y]=newVals{:};
                end
            end
            
            % Treat those that have not been matched to a ref dye, ie extra
            % target dyes
            
            w_unmatched_2=find([unmatched_sm.matched]==0);
            if ~isempty(w_unmatched_2)
                if any(sm_par(dye_id_ref).V)
                    sm_par(dye_id_ref) = get_velocity_radii(im_par, sm_par(dye_id_ref));
                    x_h=[unmatched_sm(w_unmatched_2).x];
                    y_h=[unmatched_sm(w_unmatched_2).y];
                    ds=[unmatched_sm(w_unmatched_2).diff_state];
                    sp=[unmatched_sm(w_unmatched_2).c_sp];
                    v_h = get_initial_directions(x_h, y_h, 0, sp, ds, im_par, sm_par(dye_id_ref));
                    newVals=num2cell(v_h(:,1)); [unmatched_sm(w_unmatched_2).v_x]=newVals{:};
                    newVals=num2cell(v_h(:,2)); [unmatched_sm(w_unmatched_2).v_y]=newVals{:};
                end
            end
            
            newVals = num2cell([unmatched_sm.v_x]); [sms(dye_id_target).sm(i_unmatched).v_x] = newVals{:};
            newVals = num2cell([unmatched_sm.v_y]); [sms(dye_id_target).sm(i_unmatched).v_y] = newVals{:};
            
        end
    end
else % 3D case
    for i=1:n_dyes
        % First deal with non acceptor dyes
        if ~sm_par(i).is_acceptor
            % eventually initialize velocities
            if any(sm_par(i).V)
                sm_par(i) = get_velocity_radii(im_par, sm_par(i));
                x_h=[sms(i).sm.x];
                y_h=[sms(i).sm.y];
                z_h=[sms(i).sm.z];
                ds=[sms(i).sm.diff_state];
                sp=[sms(i).sm.c_sp];
                v_h = get_initial_directions(x_h, y_h, z_h, sp, ds, im_par, sm_par(i));
                newVals=num2cell(v_h(:,1)); [sms(i).sm.v_x]=newVals{:};
                newVals=num2cell(v_h(:,2)); [sms(i).sm.v_y]=newVals{:};
                newVals=num2cell(v_h(:,3)); [sms(i).sm.v_z]=newVals{:};
            end
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
            % velocity
            [~, i_matched, i2_matched]=intersect(target_id,matched); % the ids of the matched target molecules
            newVals = num2cell([sms(dye_id_ref).sm(i2_matched).v_x]); [sms(dye_id_target).sm(i_matched).v_x] = newVals{:};
            newVals = num2cell([sms(dye_id_ref).sm(i2_matched).v_y]); [sms(dye_id_target).sm(i_matched).v_y] = newVals{:};
            newVals = num2cell([sms(dye_id_ref).sm(i2_matched).v_z]); [sms(dye_id_target).sm(i_matched).v_z] = newVals{:};
            
            % Second the sm that are not properly matched to a ref dye should
            % be assigned a diffusion state.
            [~, i_unmatched]=setdiff(target_id,matched); % the ids of the matched target molecules
            unmatched_sm=sms(dye_id_target).sm(i_unmatched);
            
            % Treat those that had been matched to a ref dye, but the ref dye
            % is not visible (fluorescent), due to maturation
            w_unmatched_1=find([unmatched_sm.matched]==-1);
            if ~isempty(w_unmatched_1)
                if any(sm_par(dye_id_ref).V)
                    sm_par(dye_id_ref) = get_velocity_radii(im_par, sm_par(dye_id_ref));
                    x_h=[unmatched_sm(w_unmatched_1).x];
                    y_h=[unmatched_sm(w_unmatched_1).y];
                    z_h=[unmatched_sm(w_unmatched_1).z];
                    ds=[unmatched_sm(w_unmatched_1).diff_state];
                    sp=[unmatched_sm(w_unmatched_1).c_sp];
                    v_h = get_initial_directions(x_h, y_h, z_h, sp, ds, im_par, sm_par(dye_id_ref));
                    newVals=num2cell(v_h(:,1)); [unmatched_sm(w_unmatched_1).v_x]=newVals{:};
                    newVals=num2cell(v_h(:,2)); [unmatched_sm(w_unmatched_1).v_y]=newVals{:};
                    newVals=num2cell(v_h(:,3)); [unmatched_sm(w_unmatched_1).v_z]=newVals{:};
                end
            end
            
            % Treat those that have not been matched to a ref dye, ie extra
            % target dyes
            
            w_unmatched_2=find([unmatched_sm.matched]==0);
            if ~isempty(w_unmatched_2)
                if any(sm_par(dye_id_ref).V)
                    sm_par(dye_id_ref) = get_velocity_radii(im_par, sm_par(dye_id_ref));
                    x_h=[unmatched_sm(w_unmatched_2).x];
                    y_h=[unmatched_sm(w_unmatched_2).y];
                    z_h=[unmatched_sm(w_unmatched_2).z];
                    ds=[unmatched_sm(w_unmatched_2).diff_state];
                    sp=[unmatched_sm(w_unmatched_2).c_sp];
                    v_h = get_initial_directions(x_h, y_h, z_h, sp, ds, im_par, sm_par(dye_id_ref));
                    newVals=num2cell(v_h(:,1)); [unmatched_sm(w_unmatched_2).v_x]=newVals{:};
                    newVals=num2cell(v_h(:,2)); [unmatched_sm(w_unmatched_2).v_y]=newVals{:};
                    newVals=num2cell(v_h(:,3)); [unmatched_sm(w_unmatched_2).v_z]=newVals{:};
                end
            end
            
            newVals = num2cell([unmatched_sm.v_x]); [sms(dye_id_target).sm(i_unmatched).v_x] = newVals{:};
            newVals = num2cell([unmatched_sm.v_y]); [sms(dye_id_target).sm(i_unmatched).v_y] = newVals{:};
            newVals = num2cell([unmatched_sm.v_z]); [sms(dye_id_target).sm(i_unmatched).v_z] = newVals{:};
            
        end
    end
end


