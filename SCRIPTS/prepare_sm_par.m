function sm_par=prepare_sm_par(n_fluorophores, sampling_rate)

sm_par(1:n_fluorophores)=struct(...
    'fluorophore_name', [], ...
    'n_mol', [], ...
    'n_mol_eff', [], ...
    'maturation_level', [], ...
    'pH_sensitivity', 0, ...
    'pH', 7.5, ...
    'spectral_data', [], ...
    'n_states', [], ...
    'state_ids', [], ...
    'initial_state', [], ...
    'initial_fluo_state', [], ...
    'fluorescent_states', [], ...
    'associated_dark_states', [], ... % In case fluorescent states are in rapid equilibrium with dark states
    'fluorescent_fraction', [], ... % In such case, fraction of the fluorescent states
    'converted_state', [], ...
    'bleached_states', [], ...
    'trans_k', [], ...
    'sampling_rate', sampling_rate, ... % sampling rate for photophysics during frametime and addtime [s-1]
    'start_sampling_rate', nan(1,int8(numel(sampling_rate)==2)+1), ... % starting sampling rate for photophysics during frametime and addtime [s-1] after optimization
    'trans_q', [], ...
    'photoactive_states', [], ...
    'quantum_yield', [], ...
    'fluorogenic', [], ...
    'fluorogenicity', [], ...
    'psf_par_ch1', [], ...
    'psf_par_ch2', [], ...
    'filter_profiles', [], ...
    'D', [], ...
    'D_ex_rates', [], ...
    'D_confined', [], ...
    'D_rate_matrix', [], ...
    'V', [], ...
    'persistence_length', [], ...
    'dispersion_selectivity', 3, ...  % power factor telling how much a sm will prefer choosing wide channels relative to narrow ones
    'margin_factor', 1e-5, ...  % Margin factor to avoid velocity vector ending up being zero when molecule close to image border
    'V_circle',  [], ...
    'V_init_dir',  [], ...
    'w_patterns', struct('w', []), ...
    'anisotropy', [], ...
    'dipole_orientation', [], ...
    'theta_fixed', [], ...
    'phi_fixed', [], ...
    'jump_allowed', 0, ... % Set to 1 if stochastic reorientation allowed
    'jump_rate', 0, ...  % Jumping rate between fixed random orientations
    'td_id', [], ... % id of tandem dye if fluorophore_pairing_on is set 
    'is_acceptor', 0, ... % if is_td=1, is dye an acceptor dye ?
    'R0_D', [], ... % Forster radius if dye is donor
    'R0_D_index', [], ... % Indices of the acceptor states able to FRET if dye is donor
    'R0_A_index', [], ... % Indices of the acceptor states able to FRET if dye is acceptor
    'processing_done', 0 ... % Will be set to 1 after dye is processed in FRET mode
    );
end