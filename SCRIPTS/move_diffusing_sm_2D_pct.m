function  sms = move_diffusing_sm_2D_pct(n_fluorophores, sms, sm_par, sm_pattern_indices, im_par, display_par)

% PURPOSE:
%	Move the molecules according to their diffusion regime in 2D
%
% INPUTS:
%   n_fluorophores: the # of dyes
%   sms: the single molecules
%	sm_par: the sm parameters
%   sm_pattern_indices: indices of virtual sample subpatterns
%	im_par: the imaging parameters
%   display_par: the displaying parameters
%
% OUTPUTS:
%   sms: the single molecules updated for: diffusion states, current position
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15
%	D.Bourgeois, September 2022, optimized for parallel computing
%	D.Bourgeois, February 2023, introduce sm_pattern_indices, now disconnected from sm_par

%Strategy:
%=> Look which molecules are in state 1, 2, 3 .... at start
%=> Look which molecules change diffusion regime at start
%=> Define new groups of molecules in state 1, 2, 3 ....
%=> Process coordinate change for each group.

if display_par.show_diff_image==1
    figure(display_par.diff_figure.Number);
end

for i=1:n_fluorophores
    %Extract the useful indices
    w_idx = find(matches(sm_par(i).sm_fn,{'x','y','x_track','y_track',...
        'sub_x','sub_y','v_x','v_y','c_sp', 'n_sp',...
        'id','bleached','activated','diff_state','n_diff_state',...
        'diff_state_trace','matched'})==1);

    sm_par(i).w_idx=w_idx; % Store this info for downstream routines

    x_idx=w_idx(1);
    y_idx=w_idx(2);
    x_track_idx=w_idx(3);
    y_track_idx=w_idx(4);
    sub_x_idx=w_idx(5);
    sub_y_idx=w_idx(6);
    id_idx=w_idx(11);
    diff_state_idx=w_idx(14);
    matched_idx=w_idx(17);

    % first look if the dye is hooked up to another dye
    if ~sm_par(i).is_acceptor
        % D: diffusion coefficients [nm^2/s]
        % C: diffusion confinement
        % K: diffusion exchange rate matrix
        D=sm_par(i).D;
        N=numel(D); % # of diffusion states
        D_state_ini=[sms(i).sm_cell{diff_state_idx,:}]; % current diffusion states of the molecules
               
        %just move those molecules that are activated and not bleached
        %in case dyes are in tandem, also consider molecules whose partner is activated and not bleached.
        w_act=sm_par(i).w_act; % Activated but not bleached molecules
        
        if im_par.fluorophore_pairing_on==1
            td_id=sm_par(i).td_id; % id of tandem molecules
            if ~isempty(td_id)
                matched=[sms(td_id).sm_cell{matched_idx,sm_par(td_id).w_act}]; % these are the activated fluorophores potentially matched to fluorophores "i"
                [~,w_act_td]=intersect([sms(i).sm_cell{id_idx,:}],matched); % this gets the index of fluorophores of type "i" that are actually matched two and activated tandem
            else
                w_act_td=[];
            end
            w_act=union(w_act, w_act_td); % these are all fluorophores to consider for moving
        end
        
        %go over all diffusion state
        for j=1:N
            w_D = find(D_state_ini==j);
            w_D_act=intersect(w_act,w_D);
            
            if ~isempty(w_D_act)
                sms(i).sm_cell(:,w_D_act) = update_diffusing_sm_2D_pct(sms(i).sm_cell(:,w_D_act), im_par, sm_pattern_indices(i), sm_par(i), j);
            end
        end
                      
        if display_par.show_diff_image==1
            scatter(im_par.m*im_par.binning-[sms(i).sm_cell{y_idx,sm_par(i).w_act}],...
                im_par.n*im_par.binning-[sms(i).sm_cell{x_idx,sm_par(i).w_act}],1,display_par.diffusion_colors(i).c(sm_par(i).w_act));
            if ~display_par.show_experiment ; pause(0.1); end % seems to be necessary for proper refresh
        end

        if im_par.fluorophore_pairing_on==1
            % we now need to move those tandem molecules that are activated
            [~,w_td_act]=intersect([sms(td_id).sm_cell{id_idx,:}],[sms(i).sm_cell{matched_idx,w_act}]); % this gets the index of fluorophores of type "i" that are actually matched to an activated tandem
            if ~isempty(w_td_act)
                pair_id=find(im_par.fluorophore_pairs(:,1)==i); % find out the pair id
                for k=1:numel(w_td_act)
                    % First identify the matched molecule
                    m=sms(td_id).sm_cell{matched_idx,w_td_act(k)};

                    % Tandem molecule will be simply moved following the
                    % attached fluorophore without consideration to
                    % confinement. This is a limitation at this stage
                    if im_par.use_diffuse_psf==1 % Initialize field sub_x sub_y
                        sms(td_id).sm_cell{sub_x_idx,w_td_act(k)} = [];
                        sms(td_id).sm_cell{sub_y_idx,w_td_act(k)} = [];
                    end
                    idx=[x_idx,y_idx,sub_x_idx,sub_y_idx];
                    sms(td_id).sm_cell(idx,w_td_act(k))=get_new_XY_td_pct(sms(i).sm_cell(idx,[sms(i).sm_cell{id_idx,:}]==m), sms(td_id).sm_cell(idx,w_td_act(k)), im_par, pair_id);
                    % Update diffusion track
                    sms(td_id).sm_cell{x_track_idx,w_td_act(k)}=vertcat(sms(td_id).sm_cell{x_track_idx,w_td_act(k)}, [sms(td_id).sm_cell{x_idx,w_td_act(k)},im_par.current_frame]);
                    sms(td_id).sm_cell{y_track_idx,w_td_act(k)}=vertcat(sms(td_id).sm_cell{y_track_idx,w_td_act(k)}, [sms(td_id).sm_cell{y_idx,w_td_act(k)},im_par.current_frame]);
                end
            end
            
            % move those unmatched and activated tandem molecules
            w_unmatched=find([sms(td_id).sm_cell{matched_idx,:}]==0);
            if ~isempty(w_unmatched)
                N=numel(td_id); % # of diffusion states
                D_state_ini=[sms(td_id).sm_cell{diff_state_idx,:}]; % current diffusion states of the molecules
                
                %just move those molecules that are unmatched and activated and not bleached
                w_act_unmatched=intersect(sm_par(td_id).w_act,w_unmatched);
                %go over all diffusion state
                for j=1:N
                    w_D = find(D_state_ini==j);
                    w_D_act=intersect(w_act_unmatched,w_D);                
                    if ~isempty(w_D_act)
                        sms(td_id).sm_cell(:,w_D_act) = update_diffusing_sm_2D_pct(sms(td_id).sm_cell(:,w_D_act), im_par, sm_pattern_indices(td_id), sm_par(td_id), j);
                    end
                end
            end
                                    
            if display_par.show_diff_image==1
                scatter(im_par.m*im_par.binning-[sms(td_id).sm_cell{y_idx,sm_par(td_id).w_act}],...
                    im_par.n*im_par.binning-[sms(td_id).sm_cell{x_idx,sm_par(td_id).w_act}],1,display_par.diffusion_colors(td_id).c(sm_par(td_id).w_act));
                if ~display_par.show_experiment ; pause(0.1); end % seems to be necessary for proper refresh
            end           
        end
    end
end