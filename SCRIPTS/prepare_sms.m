function sms=prepare_sms(n_fluorophores, sm_par, im_par)

%Define parameters for each single molecules
% sm_par(i).n_mol molecules defined for each dyes, ie incomplete maturation will be taken into account later.
sms(1:n_fluorophores)=struct('sm', []);

for i=1:n_fluorophores
    sm(1:sm_par(i).n_mol)=struct(...
        'x',0, ... % xyz position in current frame (may vary with drift) [in pixels] in high-resolution image
        'y',0, ...
        'z',0, ...
        'x_track',[], ... % list of xyz positions along the diffusion track [in pixels] in high-resolution image
        'y_track',[], ...
        'z_track',[], ...
        'sub_x',[], ... % list of xyz positions along the sub diffusion track [in pixels] if use_diffuse_psf option is used
        'sub_y',[], ...
        'sub_z',[], ...
        'v_x',[], ... % xyz velocity in current frame (may vary with drift) [in pixels] in high-resolution image if directed diffusion is used
        'v_y',[], ...
        'v_z',[], ...
        'c_sp',[], ... % # of the subpattern onto which the sm is placed along diffusion
        'n_sp',[], ... % # of the subpattern onto which the sm is placed along diffusion when/if transition occurs
        'id',0,... % id # of each molecule
        'theta',sm_par(i).theta_fixed,'phi',sm_par(i).phi_fixed,... % theta & phi polar angles. Only used if anisotropy=1
        'state', sm_par(i).initial_state, ... % Current photophysical state
        'state_id', sm_par(i).state_ids(sm_par(i).initial_state), ... % Id of current photophysical state (several states can have the same id if they interconvert rapidly and cannot be distinguished
        'state_trace', [], ... % Shows photophysical state evolution along current frame
        'fluo_trace', [], ... % Shows when the molecule is in fluorescent state along current frame
        'tot_state_trace', [0;sm_par(i).initial_state], ... % Shows photophysical state evolution along all frames until bleaching
        'tot_fluo_trace', [], ... % Shows when the molecule is in fluorescent state along all frames until bleaching
        'Ns', zeros(int8(im_par.addtime>0)+1,sm_par(i).n_states), ... % # of absorbed photons per sampling time in each starting photoactive state for addtime and frametime
        'sampling_rate', nan(1,int8(im_par.addtime>0)+1), ... % local sampling rate for addtime and frametime (might differ depending on sm position)
        'n_abs', zeros(1, sm_par(i).n_states), ... % # of absorbed photons (current frame)
        'n_em',zeros(1, sm_par(i).n_fluorescent_states), ... % # emitted photons for each fluorescent state in current frame 
        'tot_n_em',zeros(1, sm_par(i).n_fluorescent_states), ... % total # of emitted photons for each fluorescent state
        'n_phot_ch1', zeros(1, sm_par(i).n_fluorescent_states), ... % # emitted photons in ch1 for each fluorescent state (current frame)
        'n_phot_ch2', zeros(1, sm_par(i).n_fluorescent_states), ... % # emitted photons in ch2 for each fluorescent state (current frame)
        'tot_n_phot_ch1', zeros(1, sm_par(i).n_fluorescent_states), ... % total # of emitted photons in ch1 for each fluorescent state (all frames)
        'tot_n_phot_ch2', zeros(1, sm_par(i).n_fluorescent_states), ... % total # of emitted photons in ch2 for each fluorescent state(all frames)
        'n_phot_det_ch1', [], ... % total # of detected photons in ch1 (current frames)
        'n_phot_det_ch2', [], ... % total # of detected photons in ch2 (current frames)
        'tot_phot_det_ch1', 0, ... % total # of detected photons in ch1 (all frames)
        'tot_phot_det_ch2', 0, ... % total # of detected photons in ch2 (all frames)
        'frames_on_ch1', [], ... % frames where the molecule is on for channel 1 (ie, emits photons)
        'frames_on_ch2', [], ... % frames where the molecule is on for channel 2 (ie, emits photons)
        'bleached', 0, ... % 0: not bleached; 1: bleached
        'activated', 0, ... % 0: not activated; 1: activated
        'blinked', 0, ... % 0: not switched; 1: switched
        't_on', zeros(1, sm_par(i).n_fluorescent_states), ... % time the sm is in fluorescent states (current frame)
        'fr_t_on', [], ... % Fraction of times the sm is in fluorescent states (current frame) (for diffuse PSF calculations)
        't_off', zeros(1, sm_par(i).n_dark_states), ... % time the sm is in dark states (current frame)
        'diff_state', [], ... % Current Diffusion state
        'n_diff_state', [], ... % New Diffusion state when/if transition occurs
        'diff_state_trace', [;], ... % History of diffusion state
        'em_spectrum',[;], ... % Emission spectrum (current frame)
        'given_fret_photons',[], ... % Transfered # of photons for each fluorescent state when molecule is donor
        'received_fret_photons',[], ... % Transfered # of photons for each photoactive state when molecule is acceptor
        'fret_eff', 0, ... % FRET efficiency to acceptor (if donor)
        'matched',0, ... % index of partner molecule (if dyes_pair_on is set)
        'lx',[], ... % x linkage error [raster unit]
        'ly',[], ... % y linkage error [raster unit]
        'lz',[], ... % z linkage error [raster unit]
        'le_set',[] ... % 1 if linkage error is set for current frame
        );
    

    % if the molecule is not photo convertible, or if initial state is already photoconverted set all molecules to
    % activated state
    if sm_par(i).converted_state<=0
        [sm.activated]=deal(1);
    elseif sm_par(i).initial_state >= sm_par(i).converted_state
        [sm.activated]=deal(1);
    end
    newVals = num2cell(1:sm_par(i).n_mol); [sm.id] = newVals{:}; % Set id # for each molecule
    sms(i).sm=sm; clear('sm')
end