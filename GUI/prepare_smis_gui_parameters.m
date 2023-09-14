function smis_gui_parameters=prepare_smis_gui_parameters(SMIS_folder)

% SMIS vsn 2.0

smis_gui_parameters=struct(...
    'smis_title', 'My SMIS Simulation', ...
    'simul_3D', 0, ...
    'fret_on', 0, ...
    'outfiledir',pwd, ...
    'outfilename','SMIS_1', ...
    'binning',4, ...
    'raster',100, ...
    'frametime',50, ...
    'addtime',0, ...
    'n_lasers',2, ...
    'n_images',10, ...
    'n_fluorophores',1, ...
    'n_mol', [], ...
    'n_mol_eff', [], ...
    'psf_par_ch1', [], ...
    'psf_par_ch2', [], ...
    'dispersion_selectivity', 3, ...  % power factor telling how much a sm will prefer choosing wide channels relative to narrow ones
    'margin_factor', 1e-5, ...  % Margin factor to avoid velocity vector ending up being zero when molecule close to image border
    'add_diffusion',  0, ... % Set to 1 if diffusion and directed motion allowed
    'V_circle',  [], ...
    'V_init_dir',  [], ...
    'w_patterns', struct('w', []), ...
    'anisotropy', [], ...
    'dipole_orientation', [], ...
    'td_id', 0, ... % id of tandem dye if dye_pairs_on is set
    'is_acceptor', 0, ... % if is_td=1, is dye an acceptor dye ?
    'R0_D', [], ... % Forster radius if dye is donor
    'R0_D_index', [], ... % Indices of the acceptor states able to FRET if dye is donor
    'R0_A_index', [], ... % Indices of the acceptor states able to FRET if dye is acceptor
    'force_overwrite', 0, ... % Set to 1 to force overwrite files
    'processing_done', 0 ... % Will be set to 1 after dye is processed in FRET mode
    );

% Initialize default directories
if isfile(fullfile(SMIS_folder,'GUI','STARTUP_MAT_FILES','SMIS_DefaultDir.mat'))
    load(fullfile(SMIS_folder,'GUI','STARTUP_MAT_FILES','SMIS_DefaultDir.mat'),'DefaultDirs'); % Load the default directories
else
    error('SMIS Default Directory file not found !');
end

smis_gui_parameters.DefaultDirs=DefaultDirs;

% Initialize pattern files
smis_gui_parameters.in_images_dir = repmat({''},smis_gui_parameters.n_fluorophores,1);
smis_gui_parameters.in_images_names = repmat({''},smis_gui_parameters.n_fluorophores,1);

% Initialize pattern ids
smis_gui_parameters.pattern_ids(1:smis_gui_parameters.n_fluorophores)= struct(...
    'n_sp',1, ... % Number of sub patterns
    'sp_ids',[] ... % Array of sub pattern IDs
    );

% Initialize qPALM option
smis_gui_parameters.distribute_in_clusters=0; % if set to 1 n_mol/n_clusters molecules will be placed in each cluster of the image pattern
smis_gui_parameters.n_clusters=[]; % If distribute_in_clusters=1, number of clusters in the image pattern for each fluorophore

% Initialize image size
smis_gui_parameters.image_size = struct(...
    'n',[], ... % Vertical dimension
    'm',[], ... % Horizontal dimension
    'nz',[] ... % Axial dimension
    );

% Initialize drift
smis_gui_parameters.drift = struct(...
    'state',0, ... % Vertical dimension
    'x_drift',[0.2, 0, 0, 0], ... % [nm] Drift in X dimension: order 0, 1, 2 + noise [fractional]
    'y_drift',[0.2, 0, 0, 0], ... % [nm] Drift in X dimension: order 0, 1, 2 + noise [fractional]
    'z_drift',[0.2, 0, 0, 0], ... % [nm] Drift in X dimension: order 0, 1, 2 + noise [fractional]
    'rot_drift',[128, 128, 0, 0.1] ... % Rotational drift: center of rotation x0 y0 in [pixels] and rotation per frame [degrees] and noise [fractional]
    );

% Initialize objective and PSF
smis_gui_parameters.obj_and_psf = struct(...
    'obj_na',1.49, ... % Objective numerical aperture
    'obj_depth_of_focus',500, ... % Objective depth of focus
    'obj_immersion_indice', 1.515, ...  % [Water: 1.33; Oil: 1.515]
    'obj_immersion_sample', 1.33, ... % [Water: 1.33; Oil: 1.515]
    'obj_eff_from_opening_angle',1, ... % 1: Use formula from Fourkas, 2001, Opt.Letters; 0: Enter value manually
    'obj_mic_transmission', 0.8, ... % If obj_eff_from_opening_angle',1, additional correction for overall microscope transmission efficiency
    'obj_transmission_eff', 0.2, ...; % Microscope transmission efficiency, used if obj_eff_from_opening_angle=0.
    'sample_zcenter', 0, ...; % [nm] shift of sample center in Z relative to plane of focus. Set to 0 for symmetrically centered pattern. Set to -1 for plane of focus at top of sample (e.g. for TIRF mode)
    'psf_n_zslices', 32, ...; # of slices (within sample thickness)for 3D psf
    'psf_astigmatism_on', 0, ...; % Set to 1 if astigmatism for 3D SMLM is on
    'psf_astigmatism_x_save', [2.50, 1.25], ...; % Backup values for Astigmatism in X-dimension + 4rth order correction
    'psf_astigmatism_y_save', [-2.50, 1.25], ...; % Backup values for Astigmatism in Y-dimension + 4rth order correction
    'psf_astigmatism_x', [0,0], ...; % Astigmatism in X-dimension + 4rth order correction
    'psf_astigmatism_y', [0,0], ...; % Astigmatism in Y-dimension + 4rth order correction
    'psf_ch1_only', 0 ...; % Set to 1 if astigmatism to be applied only on channel 1 in two-channel experiments
    );

% Define default lasers
smis_gui_parameters.lasers(1:smis_gui_parameters.n_lasers)= struct(...
    'name','Laser', ...
    'power', 10, ...  % [mW]
    'sequence', 100*ones(1,smis_gui_parameters.n_images), ... % Percentage of laser power as a function of frame #
    'FermiPar', struct( 'NFermiFrames', smis_gui_parameters.n_images/5, ...
        'NfallofFrames', smis_gui_parameters.n_images/2), ... % Set to 'Gaussian', 'Flat'
    'RampPar', struct( 'RampStepSize', 1, ... % Step size [frames]
      'RampStart', 0, ... % Start percentage
      'RampEnd', 100), ... % End percentage
    'CustomPar', 'sin((1:n)*pi/n)', ... % Expression for custom profile
    'OnSelectedPar', '1:2', ...
    'OffSelectedPar', '1:2', ...
    'mode', 'Gaussian', ... % Set to 'Gaussian', 'Flat'
    'on_during_frametime',1, ... % Set to 1 if laser on during frametime
    'on_during_addtime', 0, ... % Set to 1 if laser on during addtime
    'beam_profile', [], ... % Laser beam profile
    'max_beam_profile', [], ... % [W/cm²] Laser maximum power density (at beam center)
    'fwhm', 50, ... % [um] size of the laser beam at sample position
    'tirf', 0, ...  % Set to 1 if 3D-TIRF mode
    'tirf_angle', 180/pi*asin(smis_gui_parameters.obj_and_psf.obj_immersion_sample/smis_gui_parameters.obj_and_psf.obj_immersion_indice), ... % [°] TIRF angle in TIRF mode; set to -1 if critical angle is to be used
    'wavelength', 500, ...  % [nm]
    'frap_mode', 0, ... % set to 1 if activation laser is restricted to a zone
    'mask_file', '', ...  % ROI File for FRAP mode
    'mask_dir', '', ...  % ROI Directory for FRAP mode
    'mask_size', [], ...  % ROI Image size for FRAP mode
    'polarization', 0, ... % [0: Circular polarization; 1: Linear polarization]
    'phi', 0 ... ; %[azimuthal angle in degrees]
    );
for i=1:smis_gui_parameters.n_lasers
    smis_gui_parameters.lasers(i).name = ['Laser ',num2str(i)];
end

% Define default fluorophore
% Load default fluorophore
load(fullfile(SMIS_folder,'GUI','STARTUP_MAT_FILES','SimplePCFP.mat'),'fluorophore');
%Load the proper fields in Fluorophorepar
f=fieldnames(fluorophore);
for i=1:numel(f)
    smis_gui_parameters.Fluorophores.(f{i})=fluorophore.(f{i});
end
%The following parameters are not intrinsic to the fluorophore so they are not defined in the fluorophore .mat file 
smis_gui_parameters.Fluorophores.N=10; % Number of fluorophores
smis_gui_parameters.Fluorophores.min_dist=4; % [nm] Minimum distance between fluorophores
smis_gui_parameters.Fluorophores.n_sp_fr=[]; % Fraction of fluorophores in each sub pattern
smis_gui_parameters.Fluorophores.Motion=struct(... %Diffusion behaviour (used if add_diffusion=1)
    'D', 0, ... % [um^2/s] Array of diffusion coefficients for  dye
    'D_rate_matrix', [], ... % Exchange rates between diffusion coefficients: array of [3 X Nx(N-1)] with N=(length(diff_coefficients)]. Array of [NxN] with N=(length(diff_coefficients)]
    'D_confined', [], ... % For each D, set to # of subpattern where diffusion is confined. Set to 0 if diffusion is confined out of subpatterns
    'D_independant_transition', 0, ... % For each D, set to 1 if a possible transition must be successful independant of where the molecule may diffuse 
    'DIT', 0, ... % Diffusion Independant transitions
    'Directed_Motion', 0, ... % Set to one if directed motion to be added
    'V', 0, ... % [um/s] Array of velocities for dye (directed diffusion)
    'persistence_length', -1); % [um] Array of persistence lengths in corresponding patterns for dye (directed diffusion only). Set to small for bulky pattern (0.01, random direction). Set to > width of pseudo-1D pattern (cytosqueletton). Set to -1 for automated choice based on chosen velocities.
    
smis_gui_parameters.Fluorophores.Orientation=struct(... % Anisotropy behavior
    'anisotropy', 0, ... % 1 means each molecule will have a fixed orientation over time
    'dipole_orientation', 1, ... %     % Used if anisotropy=1 [0: Randomly oriented , 1: Fixed orientation]
    'theta_fixed', 0, ... % Used if dipole_orientation=1 [polar angle in degrees relative to sample plane: 90°' in plane]
    'phi_fixed', 90, ... %[azimuthal angle in degrees]
    'jump_allowed', 0, ... % Set to 1 if stochastic reorientation allowed
    'jump_rate', 0);  %[Jumping rate between fixed random orientations]
  
smis_gui_parameters.Fluorophores.PSF=struct(... % Pointspread function
    'ch1', [], ... % Pointspread function of the fluorescent states in channel 1
    'ch2', []); % Pointspread function of the fluorescent states in channel 1
    

smis_gui_parameters.Fluorophore_Pairing=struct(... % fluorophore pairs to be considered for colocalization or FRET
    'state', 0, ... % Set to 1 if fluorophore pairs to be considered for colocalization or FRET
    'pairs', [], ...% for each pair, first element is id off first fluorophore, second element is id of paired fluorophore.
    'd1d2_dist', [], ...% fluorophores are covalently linked and are separated by d1d2_dist [nm]. Give one value for each pair of fluorophores 
    'd1d2_dist_sig', [], ...% fluorophores are covalently linked and are separated by d1d2_dist [nm]. Give one value for each pair of fluorophores
    'd2_constrained', []); ...% Set to 1 for tandem fluorophores that are to be positioned on their own pattern (in addition to be linked) Give one value for each pair of fluorophores

smis_gui_parameters.Fluorophores.Linkage=struct(... % fluorophore linkage property
    'set', 0, ...  % Set to 1 if linkage error allowed; only valid for FRET off
    'length', [], ... % [nm] Length of linkage error 
    'std', [], ... % [nm] Standard deviation of linkage error 
    'fixed', [], ...% Set to 1 if linkage error fixed over time; only valid for spt=0
    'pattern_control', [], ... % Set to 1 if linkage must go to defined pattern
    'sp_id', []); ...% Patterns to be reached if sm_par(k).linkage_pattern_control=1

% Define Photophysics sampling rates
smis_gui_parameters.Sampling = struct(...
    'sampling_rate',50000, ... % sampling rate for photophysics during frametime and addtime [s-1]
    'optimize_sampling_rate',1, ...  % [0/1] Set to 1 to use this option
    'min_sampling_points',5, ...  % Default minimum # of sampling points per frame
    'minimum_oversampling',10 ... % If optimize_sampling_rate==1 sampling_rate will be set to 10*max(photophysical transition rates)
    );

% Define default EMCCD parameters
smis_gui_parameters.EMCCD = struct(...
    'det_EMCCD_gain',200, ... % EMCCD gain
    'det_QE',0.9, ...  % Quantum efficiency (taken uniformly across visible spectrum). Typically a value close to 1
    'det_offset',500, ... % Detector baseline [ADU]
    'det_c',0.002, ... % Detector spurious and clock induced charge [e]
    'det_readout_noise', 74, ... % Detector rms readout noise
    'det_dynamic_range', 65535, ... % Detector dynamic range
    'det_e2p', 6, ... % Detector conversion factor: # of electrons/photon (taken uniformly across visible spectrum)
    'det_e2ADU', 45 ... % Detector gain: # of elec per ADU
    );

% Define image display parameters
smis_gui_parameters.ImageDisplay = struct(...
    'auto_display_range',1, ... % Set to 1 if automated
    'disp_lim_ch1',[0,65535], ...  % If not, give range for channel 1(min and max ADUs)
    'disp_lim_ch2',[0,65535] ...  % If not, give range for channel 1(min and max ADUs)
    );

% Define default sptPALMOptions parameters
smis_gui_parameters.sptPALMOptions = struct(...
    'move_non_activated_molecules',0, ... % Set to 1 if also the non activated molecules are moved in every frame (takes longer)
    'use_diffuse_psf',0, ...  % Set this option in case of fast diffusion so that the spot shape reflects diffusion during frametime.
    'diffuse_psf_radius',10, ... %[nm] If use_diffuse_psf=1, distance covered by diffusing molecule during frametime above which multiple PSFs will be calculated and averaged
    'ex_rates_min_oversampling',1, ... % Minimum oversampling for exchange rate calculations. Can be set << 1, eg for PAINT experiments.
    'prevent_diffusion_out_FOV',0, ... % Set to 1 to prevent diffusing molecules to escape FOV and never come back
    'show_diff_image', 0 ... % Set to 1 to see diffusion traces
    );

% Define default two-channel parameters
smis_gui_parameters.TwoChannel = struct(...
    'state',0, ... % Set to 1 to enable second channel
    'defocus', 1, ...  % 1: No defocus; psf width will be multiplied by two_channel_defocus
    'deform', struct(... % Second channel deformation
    'X_shift', 0, ... % X-shift[raster]
    'Y_shift', 0, ... % Y-shift[raster]
    'X_Stretch', 1, ... % X-stretch [fraction over whole image]
    'Y_Stretch', 1, ... % Y-stretch [fraction over whole image]
    'Rotation', 0, ... % rotation [deg] from image center
    'Noise', 0), ... % noise [in %]
    'single_CCD',0, ... % Set to 1 if the 2 channels are recorded on a single CCD
    'em_filter_1', '[500 40 1]', ... % Emission filter for channel 1
    'use_exp_ch1_filter', 0, ... % Set to 1 if commercial filter is used
    'exp_ch1_filter', [], ... % Commercial channel 1 filter
    'em_filter_2', '[600 40 1]', ... % Emission filter for channel 2
    'use_exp_ch2_filter', 0, ... % Set to 1 if commercial filter is used
    'exp_ch2_filter', [], ... % Commercial channel 2 filter
    'add_dichroic', 0, ... % Set to 1 if Dichroic filter is used
    'dichroic_filter', '[600 10 0.05 0.95]', ... % Dichroic filter
    'use_exp_dichroic', 0, ... % Set to 1 if commercial dichroic filter is used
    'exp_dichroic_filter', [], ... % Commercial dichroic filter
    'filter_path', [] ... % Path for commercial filters
    );

default_bg_laser_sensitivity=(1/smis_gui_parameters.n_lasers)*ones(1,smis_gui_parameters.n_lasers);
smis_gui_parameters.BG = struct(...
    'bg_ch1', 1, ... % intensity [photons/100x100nm^2/frame] for channel 1, Poisson noise is assumed
    'bg_ch2', 1, ... % intensity [photons/100x100nm^2/frame] for channel 2, Poisson noise is assumed
    'adjust_bg', 0,... % [1:Yes/0:No]
    'bg_laser_sensitivity_ch1',default_bg_laser_sensitivity, ... % if adjust_bg=1, fractional level of laser sensitivity (total must be 1)
    'bg_laser_sensitivity_ch2',default_bg_laser_sensitivity, ... % if adjust_bg=1, fractional level of laser sensitivity (total must be 1)
    'bg_decay_rate',0.1, ... % [(uJ/100x100nm^2)-1] if adjust_bg=1, decay rate of background fluorescence
    'add_textured_bg',0, ... % [1:Yes/0:No]
    'textured_bg_pattern_dir','', ...
    'textured_bg_pattern_file', '', ...
    'textured_bg_ch1', 50, ... % intensity [photons/100x100nm^2/frame] for channel 1, Poisson noise is assumed
    'textured_bg_ch2', 50, ... % intensity [photons/100x100nm^2/frame] for channel 2, Poisson noise is assumed
    'adjust_textured_bg', 0, ... % [1:Yes/0:No]
    'textured_bg_laser_sensitivity_ch1',1/smis_gui_parameters.n_lasers*ones(1,smis_gui_parameters.n_lasers), ... % if adjust_textured_bg=1, level of laser sensitivity between 0 and 1
    'textured_bg_laser_sensitivity_ch2',1/smis_gui_parameters.n_lasers*ones(1,smis_gui_parameters.n_lasers), ... % if adjust_textured_bg=1, level of laser sensitivity between 0 and 1
    'textured_bg_decay_rate', 0 ... % [(uJ/100x100nm^2)-1] if adjust_textured_bg=1, decay rate of textured background fluorescence
    );

% Define general parameters
smis_gui_parameters.General = struct(...
    'stochastic_spectrum_off',1, ... % To avoid stochastic spectral distribution of emitted photons and save a lot of time turn this option on ([1:Yes, 0:No])
    'apply_poisson_stat_for_lasers', 0, ...  % Set to 1 to treat # of absorbed laser photons with Poisson statistics
    'use_parallel_computing', 0, ... % Set to 1 to use parallel computing
    'min_n_mol_pct', 1000 ... % Minimum # of molecule for parallel computing to be useful
    );

% Define displaying options
smis_gui_parameters.show_experiment=0;
smis_gui_parameters.show_lasers=0; % set to 1 if laser profiles and intensity to be shown
smis_gui_parameters.show_dye_spectra=0; % set to 1 if spectra from dyes to be shown
smis_gui_parameters.show_psf=0; % set to 1 if psf to be shown
smis_gui_parameters.show_photophysics=0; % set to 1 if photophysics to be shown at the end of simulation
smis_gui_parameters.show_data_aquisition=0; % set to 1 if data aquisition to be shown
smis_gui_parameters.keep_current_figure=1; % set to 1 if data aquisition to be shown on current figure
smis_gui_parameters.lasers_figure_number=[]; % Handle to the figure number displaying lasers
smis_gui_parameters.spectra_figure_number=[]; % Handle to the figure number displaying spectra
smis_gui_parameters.emccd_figure_number=[]; % Handle to the figure number displaying lasers
smis_gui_parameters.photophysics_figure_number=[]; % Handle to the figure number displaying photophysics summary
smis_gui_parameters.psf_figure_number=[]; % Handle to the figure number displaying spectra
end

