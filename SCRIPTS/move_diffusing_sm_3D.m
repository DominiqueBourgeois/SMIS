function  sms = move_diffusing_sm_3D(n_fluorophores, sms, im_par, sm_par, display_par)

% PURPOSE:
%	Move the molecules according to their diffusion regime in 3D
%
% INPUTS:
%   n_fluorophores: the # of fluorophores
%   sms: the single molecules (with coordinates on high-resolution image)
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   sms: the single molecules updated for: diffusion states, current position
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15

%Strategy:
%=> Look which molecules are in state 1, 2, 3 .... at start
%=> Look which molecules change diffusion regime at start
%=> Define new groups of molecules in state 1, 2, 3 ....
%=> Process coordinate change for each group.

if display_par.show_diff_image==1
    figure(display_par.diff_figure.Number);
end

for i=1:n_fluorophores
    % first look if the dye is hooked up to another dye
    if ~sm_par(i).is_acceptor
        % D: diffusion coefficients [nm^2/s]
        % C: diffusion confinement
        % K: diffusion exchange rate matrix
        D=sm_par(i).D;
        N=numel(D); % # of diffusion states
        D_state_ini=[sms(i).sm.diff_state]; % current diffusion states of the molecules
        
        %just move those molecules that are activated and not bleached
        %in case dyes are in tandem, also consider molecules whose partner is activated and not bleached.
        w_act=sm_par(i).w_act;

        if im_par.fluorophore_pairing_on==1
            td_id=sm_par(i).td_id; % id of tandem molecules
            if ~isempty(td_id)
                matched=[sms(td_id).sm(sm_par(td_id).w_act).matched]; % these are the activated fluorophores potentially matched to fluorophores "i"
                [~,w_act_td]=intersect([sms(i).sm.id],matched); % this gets the index of fluorophores of type "i" that are actually matched two and activated tandem
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
                sms(i).sm(w_D_act) = update_diffusing_sm_3D(sms(i).sm(w_D_act), im_par, sm_par(i), j);
            end
        end
                
        if display_par.show_diff_image==1
            scatter3(im_par.m*im_par.binning-[sms(i).sm(sm_par(i).w_act).y],...
                im_par.n*im_par.binning-[sms(i).sm(sm_par(i).w_act).x],...
                im_par.nz*im_par.binning-[sms(i).sm(sm_par(i).w_act).z], 1,display_par.diffusion_colors(i).c(sm_par(i).w_act));
            if ~display_par.show_experiment ; pause(0.1); end % seems to be necessary for proper refresh
        end
        
        if im_par.fluorophore_pairing_on==1
            % we now need to move those tandem molecules that are activated
            [~,w_td_act]=intersect([sms(td_id).sm.id],[sms(i).sm(w_act).matched]); % this gets the index of fluorophores of type "i" that are actually matched two and activated tandem
            if ~isempty(w_td_act)
                pair_id=find(im_par.fluorophore_pairs(:,1)==i); % find out the pair id
                for k=1:numel(w_td_act)
                    % First identify the matched molecule
                    m=sms(td_id).sm(w_td_act(k)).matched;
                    
                    % Tandem molecule will be simply moved following the
                    % attached fluorophore without consideration to
                    % confinement. This is a limitation at this stage
                    if im_par.use_diffuse_psf==1 % Initialize field sub_x sub_y sub_z
                        sms(td_id).sm(w_td_act(k)).sub_x = [];
                        sms(td_id).sm(w_td_act(k)).sub_y = [];
                        sms(td_id).sm(w_td_act(k)).sub_z = [];
                    end
                    sms(td_id).sm(w_td_act(k))=get_new_XYZ_td(sms(i).sm([sms(i).sm.id]==m), sms(td_id).sm(w_td_act(k)), im_par, pair_id);
                    % Update diffusion track
                    sms(td_id).sm(w_td_act(k)).x_track=vertcat(sms(td_id).sm(w_td_act(k)).x_track, [sms(td_id).sm(w_td_act(k)).x,im_par.current_frame]);
                    sms(td_id).sm(w_td_act(k)).y_track=vertcat(sms(td_id).sm(w_td_act(k)).y_track, [sms(td_id).sm(w_td_act(k)).y,im_par.current_frame]);
                    sms(td_id).sm(w_td_act(k)).z_track=vertcat(sms(td_id).sm(w_td_act(k)).z_track, [sms(td_id).sm(w_td_act(k)).z,im_par.current_frame]);
                end
            end
            
            % move those unmatched and activated tandem molecules
            w_unmatched=find([sms(td_id).sm.matched]==0);
            if ~isempty(w_unmatched)
                N=numel(td_id); % # of diffusion states
                D_state_ini=[sms(td_id).sm.diff_state]; % current diffusion states of the molecules
                
                %just move those molecules that are unmatched and activated and not bleached
                w_act_unmatched=intersect(sm_par(td_id).w_act,w_unmatched);
                %go over all diffusion state
                for j=1:N
                    w_D = find(D_state_ini==j);
                    w_D_act=intersect(w_act_unmatched,w_D);
                    if ~isempty(w_D_act)
                        sms(td_id).sm(w_D_act) = update_diffusing_sm_3D(sms(td_id).sm(w_D_act), im_par, sm_par(td_id), j);
                    end
                end
            end
            
            if display_par.show_diff_image==1
                scatter3(im_par.m*im_par.binning-[sms(td_id).sm(sm_par(td_id).w_act).y],...
                    im_par.n*im_par.binning-[sms(td_id).sm(sm_par(td_id).w_act).x],...
                    im_par.nz*im_par.binning-[sms(td_id).sm(sm_par(td_id).w_act).z], 1,display_par.diffusion_colors(td_id).c(sm_par(td_id).w_act));
                if ~display_par.show_experiment ; pause(0.1); end % seems to be necessary for proper refresh
            end
        end
    end
end
