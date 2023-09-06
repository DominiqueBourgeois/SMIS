% Single Molecule Imaging Simulator (SMIS) 
% Ensemble case
% Dominique Bourgeois
% December, 2021

function smis_ok=smis_vsn2_2_ensemble_main(app, smis_par)

% clear all;
% dbstop if error % Stop debugger inside called functions
vsn_no='2.2'; % Version number

if ~isdeployed % Take this line off for standalon version
    %Path where the simulation software is loacted
    addpath(genpath(smis_par.DefaultDirs.SMIS_GUI_Dir(1:end-4)));
end


% Number of fluorophores
n_fluorophores=smis_par.n_fluorophores; % # of fluorophores

% simulation in 3-D ?
simul_3D=smis_par.simul_3D; % [Yes=1, No=0] Set to 1 for simulations in 3D

%For each fluorophore read a pattern image: the image should be a .tif in 2D and a 3D-kernel in 3D
in_images_dir = smis_par.in_images_dir;
in_images_names = smis_par.in_images_names;

%number of clusters in cluster pattern
distribute_in_clusters=smis_par.distribute_in_clusters; % if set to 1 n_mol/n_clusters molecules will be placed in each cluster of the image pattern
n_clusters=smis_par.n_clusters; % If distribute_in_clusters=1, number of clusters in the image pattern for each fluorophore

% Binning factor: the input (or created) image is a high-resolution image and is "binning"^2
% larger in size than the (virtual) detected image and its pixel size is "binning"
% smaller in each dimension. In 3D the Z-pixel size in the high-resolution
% kernel is also set to: raster/binning
binning=smis_par.binning;

%2D raster size [nm] (of the detected image).
raster=smis_par.raster; % [nm]

%Frametime [ms] (camera on time).
frametime=smis_par.frametime; % [ms]

%addtime [ms] (additionnal time before camera on time in each cycle).
addtime=smis_par.addtime; % [ms]

%sampling rate [Hz]: fluorescence traces will the computed at this rate at
%a minimum. Trace resolution, but also computing time will scale with this rate.
%the sampling rate might be increased if necessary by the software if
%required by the fluorescent state transition matrix for the various fluorophores.
%if 2 numbers are provided, the first will be used during frametime and the
%second during addtime.
sampling_rate=smis_par.Sampling.sampling_rate; % 1e+04; % [Hz]
if addtime>0 
    sampling_rate=[sampling_rate,sampling_rate]; % Use the same sampling rate for addtime to start with
end

%It is also possible to optimize sampling rate so as to minimize computing time
optimize_sampling_rate=smis_par.Sampling.optimize_sampling_rate; % [0/1] Set to 1 to use this option
minimum_oversampling=smis_par.Sampling.minimum_oversampling; % If optimize_sampling_rate==1 sampling_rate will be set to 10*max(photophysical transition rates)
min_sampling_points=smis_par.Sampling.min_sampling_points; % Minimum # of sampling points per frametime or addtime

%is it a two-channel experiment (1: Yes, 0: No)
%Yes means dual-channel (possibly only one fluorophore)
two_channel=smis_par.TwoChannel.state;


two_channel_defocus=smis_par.TwoChannel.defocus; % (1: No defocus; psf width will be multiplied by two_channel_defocus)

two_channel_deform= ...
    [smis_par.TwoChannel.deform.X_shift, smis_par.TwoChannel.deform.X_shift, ... % X-shift[raster], Y-shift
    smis_par.TwoChannel.deform.X_Stretch, smis_par.TwoChannel.deform.Y_Stretch, ... % X-stretch [fraction over whole image], Y-stretch
    smis_par.TwoChannel.deform.Rotation, smis_par.TwoChannel.deform.Noise]; % rotation [deg] from image center, noise [in %])

%If two_channel=1, are the 2 channels recorded on a single CCD
single_CCD=smis_par.TwoChannel.single_CCD; % (1: Yes, 0: No)

%define channels by filters; Multiband filters must be with bands of increasing values
ch1_filter = eval(smis_par.TwoChannel.em_filter_1); % [Low wavelength [nm], High wavelength [nm], transmission efficiency] or [Central wavelength [nm], width [nm], transmission efficiency]
ch2_filter = eval(smis_par.TwoChannel.em_filter_2); % [Low wavelength [nm], High wavelength [nm], transmission efficiency] or [Central wavelength [nm], width [nm], transmission efficiency]
use_exp_ch1_filter=smis_par.TwoChannel.use_exp_ch1_filter; % Set to 1 if experimental (commercial) filter profile is applied for channel 1
use_exp_ch2_filter=smis_par.TwoChannel.use_exp_ch2_filter; % Set to 1 if experimental (commercial) filter profile is applied for channel 2
exp_ch1_filter = smis_par.TwoChannel.exp_ch1_filter; % Experimental profile from commercial filter, chanel 1
exp_ch2_filter = smis_par.TwoChannel.exp_ch2_filter; % Experimental profile from commercial filter, chanel 2
%dichroic filter
add_dichroic=smis_par.TwoChannel.add_dichroic; % smis_par.TwoChannel.add_dichroic;
use_exp_dichroic=smis_par.TwoChannel.use_exp_dichroic; % smis_par.TwoChannel.use_exp_dichroic;
dichroic_filter = eval(smis_par.TwoChannel.dichroic_filter); % eval(smis_par.TwoChannel.dichroic_filter); % [Cutoff wavelength [nm], cutoff width [nm], transmission efficiency low pass in ch2, transmission efficiency low pass in ch1]
exp_dichroic_filter=smis_par.TwoChannel.exp_dichroic_filter; % Experimental profile from dichroic filter

%handle point-spread-function (PSF)
psf_mode='Simple_Gaussian'; % Can be obtained from 'Simple_Gaussian', 'Wavefront Toolbox', 'EPFL_psf' or 'scalar diffraction-limited model'
% for now only 'Simple_Gaussian' is implemented
% Only for 3D; for now used with 'Simple_Gaussian' psf_mode
% See Huang et al Biomed Opt Express, 2015, 6, 902

psf_astigmatism_x=smis_par.obj_and_psf.psf_astigmatism_x; % Astigmatism in X-dimension + 4rth order correction
psf_astigmatism_y=smis_par.obj_and_psf.psf_astigmatism_y; % Astigmatism in Y-dimension + 4rth order correction
psf_astigmatism_ch1_only=smis_par.obj_and_psf.psf_ch1_only; % Set to 1 if astigmatic PSF to be used only in ch1
psf_n_zslices=smis_par.obj_and_psf.psf_n_zslices; % # of slices (within sample thickness)for 3D psf
sample_zcenter=smis_par.obj_and_psf.sample_zcenter; % [nm] shift of sample center in Z relative to plane of focus. Set to 0 for symmetrically centered pattern. Set to -1 for plane of focus at top of sample (e.g. for TIRF mode)
depth_of_focus=smis_par.obj_and_psf.obj_depth_of_focus; % [nm] Depth of focus of objective/microscope setup. Set to -1 for automated calculation [to be implemented]

n_images=smis_par.n_images; % if yes, # of images
outfilename=smis_par.outfilename;
outfiledir=smis_par.outfiledir;
force_overwrite=smis_par.force_overwrite;

%Add drift ? ([1:Yes, 0:No])
add_drift=smis_par.drift.state;

%To avoid stochastic spectral distribution of emitted photons and save
%a lot of time turn this option on ([1:Yes, 0:No])
stochastic_spectrum_off=smis_par.General.stochastic_spectrum_off;

%Define lasers
n_lasers=smis_par.n_lasers; % Number of lasers
lasers=prepare_lasers(n_lasers); % Define the lasers
apply_poisson_stat_for_lasers=smis_par.General.apply_poisson_stat_for_lasers; % Set to 1 to treat # of absorbed laser photons with Poisson statistics

for k=1:n_lasers
    %Readout laser parameters
    lasers(k).name=smis_par.lasers(k).name;
    lasers(k).power = smis_par.lasers(k).power; % [mW]
    lasers(k).sequence =smis_par.lasers(k).sequence; % Constant. Percentage of laser power as a function of frame #
    lasers(k).on_during_frametime = smis_par.lasers(k).on_during_frametime; % Set to 1 if laser on during frametime
    lasers(k).on_during_addtime = smis_par.lasers(k).on_during_addtime; % Set to 1 if laser on during addtime
    lasers(k).mode=smis_par.lasers(k).mode; % Set to 'Gaussian', 'Flat'
    lasers(k).fwhm = smis_par.lasers(k).fwhm; % [um] fwhm of the laser beam at focus; valid for both 'Gaussian' and 'Flat'
    lasers(k).tirf = smis_par.lasers(k).tirf; % Set to 1 if 3D-TIRF or 3D-HILO mode
    lasers(k).tirf_angle = smis_par.lasers(k).tirf_angle; % [°] TIRF angle in TIRF mode; set to -1 if critical angle is to be used ; if angle is below critical angle, HILO mode will be used
    lasers(k).wavelength = smis_par.lasers(k).wavelength; % [nm]
    %If anisotropy=1, check laser polarization
    %Polarization is assumed perpendicular to beam axis
    lasers(k).polarization=smis_par.lasers(k).polarization; % [0: Circular polarization; 1: Linear polarization]
    %If laser_polarization=1;
    lasers(k).phi=smis_par.lasers(k).phi; %[azimuthal angle in degrees]
    lasers(k).frap_mode = smis_par.lasers(k).frap_mode; % set to 1 if activation laser is restricted to a zone
    lasers(k).mask_file = smis_par.lasers(k).mask_file;
    lasers(k).mask_dir=smis_par.lasers(k).mask_dir;
end


%Add diffusion ? ([1:Yes, 0:No])
add_diffusion=smis_par.add_diffusion;
move_non_activated_molecules=smis_par.sptPALMOptions.move_non_activated_molecules; % Set to 1 if also the non activated molecules are moved in every frame (takes longer)
use_diffuse_psf=smis_par.sptPALMOptions.use_diffuse_psf; % Set this option in case of fast diffusion so that the spot shape reflects diffusion during frametime.
diffuse_psf_radius=smis_par.sptPALMOptions.diffuse_psf_radius; %[nm] If use_diffuse_psf=1, distance covered by diffusing molecule during frametime above which multiple PSFs will be calculated and averaged
prevent_diffusion_out_FOV=smis_par.sptPALMOptions.prevent_diffusion_out_FOV; % Set to 1 to prevent diffusing molecules to escape FOV and never come back
show_diff_image = smis_par.sptPALMOptions.show_diff_image; % Set to 1 to see diffusion traces

% Definition of fluorophores
sm_par=prepare_sm_par(n_fluorophores, sampling_rate);

% Definition of populations
ens(1:n_fluorophores)=struct;

% Go through all fluorophores
for k=1:n_fluorophores
    %fluorophore name
    sm_par(k).fluorophore_name = smis_par.Fluorophores(k).name;
    % Number of molecules
    sm_par(k).n_mol=smis_par.Fluorophores(k).N;
    
    % The molecules may be distributed in different
    % subpatterns, depending on the enterered pattern images. Give the
    % subpattern ids, ie values of segmented subpatterns.
    % 1 generally corresponds to the first subpattern, 2 to the second, and 0 must correspond to elsewhere in the image
    sm_par(k).n_sp_id = smis_par.pattern_ids(k).sp_ids';
    %In that case, the fluorophores are distributed between the subpatterns
    sm_par(k).n_sp_fr = smis_par.Fluorophores(k).n_sp_fr;
    
    % Maturation level
    sm_par(k).maturation_level=smis_par.Fluorophores(k).maturation_level; % Set to < 1 for incomplete maturation level (n_mol*maturation_level molecules will be used)
    
    %Nearest neighbour distance between molecules [nm]
    sm_par(k).near_dist = smis_par.Fluorophores(k).min_dist; % [nm] 2 molecules cannot end up being distant by less than near_dist
    
    % pH
    sm_par(k).pH_sensitivity=smis_par.Fluorophores(k).Environment.pH_sensitivity; % Set to 1 if fluorophore is pH sensitive. If set to 0 the fluorescent states must not be in equilibrium with non-fluorescent states
    sm_par(k).pH=smis_par.Fluorophores(k).Environment.pH; % pH of the medium for the fluorophore
    
    sm_par(k).state_names=smis_par.Fluorophores(k).Photophysics.state_names; % Names of the photophysical states
    sm_par(k).n_states=smis_par.Fluorophores(k).Photophysics.n_states;   % Total number of photophysical states
    sm_par(k).state_ids=smis_par.Fluorophores(k).Photophysics.state_ids; % ids of all states for the sm: BE CAREFUL: give the same id for states that are in rapid equilibrium
    sm_par(k).initial_state=sm_par(k).state_ids(smis_par.Fluorophores(k).Photophysics.initial_state); % # of the initial state for the sm: Take care of the states in rapid equilibrium
    sm_par(k).fluorescent_states=smis_par.Fluorophores(k).Photophysics.fluorescent_states; % # of the fluorescent state for the sm
    sm_par(k).converted_state=smis_par.Fluorophores(k).Photophysics.converted_state; % # of the activated state for the sm
    sm_par(k).bleached_states=smis_par.Fluorophores(k).Photophysics.bleached_states; % # of the bleached state for the sm
    sm_par(k).photoactive_states=smis_par.Fluorophores(k).Photophysics.photoactive_states; % # of the photoactive state for the sm
    
    %Spectral data
    sm_par(k).spectral_data.n_fluorescent_states=smis_par.Fluorophores(k).Photophysics.n_fluorescent_states;
    sm_par(k).spectral_data.n_photoactive_dark_states=smis_par.Fluorophores(k).Photophysics.n_photoactive_dark_states;
    sm_par(k).spectral_data.em_spectra=smis_par.Fluorophores(k).Photophysics.em_spectra;
    sm_par(k).spectral_data.exc_spectra=smis_par.Fluorophores(k).Photophysics.exc_spectra;
    sm_par(k).spectral_data.dark_spectra=smis_par.Fluorophores(k).Photophysics.dark_spectra;
    
    
    %Thermal transitions (give rates in [s-1], rates should not exceed 100kHz at most)
    sm_par(k).trans_k=smis_par.Fluorophores(k).Photophysics.trans_k; %  thermally activated transition rate  matrice
    %Photoactivated transitions (give quantum yields)
    sm_par(k).trans_q=smis_par.Fluorophores(k).Photophysics.trans_q; %  photoactivated transition QY  matrice
    
    % Fluorescence quantum yield:  give one value for each fluorescent state
    sm_par(k).quantum_yield =smis_par.Fluorophores(k).Photophysics.quantum_yield;
    
    %pKa and Hill coefficient
    sm_par(k).pKa =smis_par.Fluorophores(k).Environment.pKa;
    sm_par(k).Hill =smis_par.Fluorophores(k).Environment.Hill;
    
    %handle anisotropy (1: Yes (typically fixed cell); 0: No (typically live cell))
    sm_par(k).anisotropy=smis_par.Fluorophores(k).Orientation.anisotropy; % 1 means each molecule will have a fixed orientation over time
    %If anisotropy=1, use the following parameters
    sm_par(k).dipole_orientation=smis_par.Fluorophores(k).Orientation.dipole_orientation; % [0: Randomly _oriented dipole; 1: Fixed orientation]
    %If dipole_orientation=1;
    sm_par(k).theta_fixed=smis_par.Fluorophores(k).Orientation.theta_fixed; %[polar angle in degrees relative to sample plane: 90° = in plane]
    sm_par(k).phi_fixed=smis_par.Fluorophores(k).Orientation.phi_fixed; %[azimuthal angle in degrees]
        
    %FOR ENSEMBLE SIMULATIONS
    sm_par(k).n_subpop=1; % Define just 1 sub population for now TO BE DEVELOPPED: SPLIT IN BINS 
    
    ens(k).state_names=smis_par.Fluorophores(k).Photophysics.state_names;
    ens(k).n_states=smis_par.Fluorophores(k).Photophysics.n_states;
    ens(k).sp(1:sm_par(k).n_subpop)=struct;
    for i=1:sm_par(k).n_subpop
        ens(k).sp(i).p=zeros(1,ens(k).n_states); % Current population of the states over time
        ens(k).sp(i).p(1,sm_par(k).initial_state)=1; %Set up initial state       
        ens(k).sp(i).det_p=zeros(ens(k).n_states,smis_par.n_images); % Detected subpopulations of the states over collected frames
        ens(k).sp(i).S_ch1=zeros(smis_par.Fluorophores(k).Photophysics.n_fluorescent_states,smis_par.n_images); % Detected populations of the states over collected frames
        if smis_par.TwoChannel.state==1
            ens(k).sp(i).S_ch2=zeros(smis_par.Fluorophores(k).Photophysics.n_fluorescent_states,smis_par.n_images); % Detected populations of the states over collected frames
        end
        ens(k).sp(i).x=[]; % x position
        ens(k).sp(i).y=[]; % y position
        ens(k).sp(i).z=[]; % z position
        ens(k).sp(i).sampling_rate=nan(int8(smis_par.addtime>0)+1,1); % local sampling rate for addtime and frametime (might differ depending on sm position)
        ens(k).sp(i).Ns=zeros(int8(smis_par.addtime>0)+1,ens(k).n_states); % # of absorbed photons per sampling time in each starting photoactive state for addtime and frametime
    end
     ens(k).det_p=zeros(size(ens(k).sp(1).det_p)); % Detected total populations of the states over collected frames
     ens(k).S_ch1=zeros(size(ens(k).sp(1).S_ch1));
     if smis_par.TwoChannel.state==1
         ens(k).S_ch2=zeros(size(ens(k).sp(1).S_ch2));
     end
     
end

%Tandem fluorophores
fluorophore_pairing_on=smis_par.Fluorophore_Pairing.state; % Set to 1 if fluorophore pairs to be considered for colocalization or FRET
fluorophore_pairs=smis_par.Fluorophore_Pairing.pairs; % for each pair, first element is sm_par id, second element is paired sm_par id.
%If FRET is turned on, first element is DONOR and second element is ACCEPTOR
%Second element fluorophores will be positionned in same patterns as first element
%fluorophores, except additionnal second element fluorophores that my be unpaired
%Average distance between pairs of fluorophores nm]
d1d2_dist=smis_par.Fluorophore_Pairing.d1d2_dist ; % fluorophores are covalently linked and are separated by d1d2_dist [nm]. Give one value for each pair of fluorophores
%fluctuation for distance between pairs of fluorophores [nm]
d1d2_dist_sig=smis_par.Fluorophore_Pairing.d1d2_dist_sig;  % distance between fluorophores fluctuate from one pair to another by this amount. Give one value for each pair of fluorophores
d2_constrained=smis_par.Fluorophore_Pairing.d2_constrained; % Set to 1 for tandem fluorophores that are to be positioned on their own pattern (in addition to be linked) Give one value for each pair of fluorophores

%Turn FRET between the tandems on or off ([1:On, 0:Off])
fret_on=smis_par.fret_on;


% Objective parameters
% PSF will be assumed to be symmetrical Gaussian with fwhm = 1.22*lambda/2NA
% with lambda =  laser wavelength
obj_na = smis_par.obj_and_psf.obj_na;
% Indice of refraction of immersion medium
obj_immersion_indice=smis_par.obj_and_psf.obj_immersion_indice; % [Water: 1.33; Oil: 1.515]
% Indice of refraction of sample
obj_immersion_sample=smis_par.obj_and_psf.obj_immersion_sample; % [Water: 1.33; Oil: 1.515]

% detection efficiency [fraction of emitted photons detected by objective and optical path]
obj_eff_from_opening_angle=smis_par.obj_and_psf.obj_eff_from_opening_angle; % 1: Use formula from Fourkas, 2001, Opt.Letters; 0: Enter value manually
obj_mic_transmission=smis_par.obj_and_psf.obj_mic_transmission; % If obj_eff_from_opening_angle=1, additional correction for overall microscope transmission efficiency
obj_transmission_eff = smis_par.obj_and_psf.obj_transmission_eff; % Microscope transmission efficiency, used if obj_eff_from_opening_angle=0.

% Show experiment in various figures
show_experiment=smis_par.show_experiment;
keep_current_figure=smis_par.keep_current_figure; % set to 1 if data aquisition to be shown on current figure

%Turn on debug mode [Yes=1; No=0]
debug=0;

%Variable to control stopping of simulation from SMIS GUI
smis_interrupt=0; % Use that if simulation is interrupted !
exec_time=0.1; % Initialize execution time for one image

%   *********** END OF INPUT ***********

%% Start
%Display starting message and image size
clc; % clear command window
% Display title
disp('****************************************************');
disp(['Running SMIS simulation: ',smis_par.smis_title, '  ...']);
disp(['********           vsn:',vsn_no,'        **********']);
disp('******** ENSEMBLE SIMULATION **********');
disp('****************************************************');

if (show_experiment==1 || show_diff_image) && keep_current_figure==0
    close all;
end

clear('a_h_all','a_h_all_sm','sms','im_par','photophysical_status','photophysical_status_2','photophysical_status_td','photophysical_status_2_td')

%% Fixes and restrictions

if fret_on==1
    SMISMessage='FRET not available in ensemble simulations !';
    disp(SMISMessage);
    warndlg(SMISMessage,'Warning')
    smis_ok=0; return
end

if fluorophore_pairing_on==1 && max(fluorophore_pairs(:))>n_fluorophores
    SMISMessage=['fluorophore_pairs is not properly set for only ', num2str(n_fluorophores), ' fluorophores !'];
    disp(SMISMessage);
    warndlg(SMISMessage,'Warning')
    smis_ok=0; return
end

if smis_par.Sampling.optimize_sampling_rate==1
    SMISMessage='Optimization of sampling rate not implemented in ensemble simulations !';
    disp(SMISMessage);
    warndlg(SMISMessage,'Warning')
    smis_ok=0; return
end


if size(sm_par,2)>n_fluorophores % In case too many fluorophores defined
    sm_par(n_fluorophores+1:end)=[];
end

%Check for pH sensitivity
for i=1:n_fluorophores
    if sm_par(i).pH_sensitivity==0 % There must not be dark states in rapid equilibrium with fluorescent states
        for j=1:numel(sm_par(i).fluorescent_states)
            % Is there a state in equilibrium ?
            if sum(sm_par(i).state_ids==sm_par(i).fluorescent_states(j))>1
                SMISMessage1=['Fluorophore ', num2str(i), ': Fluorescent state ', num2str(sm_par(i).fluorescent_states(j)), ' has a dark state in rapid equilibrium ! '];
                SMISMessage2='This is not possible when pH_sensitivity is set to 0 ! Fix fluorophore definition ';
                SMISMessage=[SMISMessage1 newline SMISMessage2];
                disp(SMISMessage);
                warndlg(SMISMessage,'Warning')
                smis_ok=0; return
            end
        end
    else
        for j=1:numel(sm_par(i).fluorescent_states)
            % Is there a state in equilibrium ?
            if sum(sm_par(i).state_ids==sm_par(i).fluorescent_states(j))==1
                SMISMessage1=['Fluorophore ', num2str(i), ': Fluorescent state ', num2str(sm_par(i).fluorescent_states(j)), ' has no dark state in rapid equilibrium ! '];
                SMISMessage2='This is not possible when pH_sensitivity is set to 1 ! Fix fluorophore definition ';
                SMISMessage=[SMISMessage1 newline SMISMessage2];
                disp(SMISMessage);
                warndlg(SMISMessage,'Warning')
                smis_ok=0; return
            end
        end
    end
    
    if sm_par(i).fluorogenic==1
        SMISMessage=['Fluorogenicity not handled in ensemble simulations ! Fluorophore # ', num2str(i)];
        disp(SMISMessage);
        warndlg(SMISMessage,'Warning')
        smis_ok=0; return
    end
    
    if sm_par(i).anisotropy
        SMISMessage=['Anisotropy not handled in ensemble simulations ! Fluorophore # ', num2str(i)];
        disp(SMISMessage);
        warndlg(SMISMessage,'Warning')
        smis_ok=0; return
    end
end

if any([sm_par.maturation_level]>1)
    SMISMessage='Maturation level cannot exceed 1!';
    disp(SMISMessage);
    warndlg(SMISMessage,'Warning')
    smis_ok=0; return
end

if simul_3D==1 && any([lasers.tirf]==1)
    SMISMessage='TIRF mode not available in ensemble simulations';
    disp(SMISMessage);
    warndlg(SMISMessage,'Warning')
    smis_ok=0; return
end

if add_drift==1
    SMISMessage='Drift not available in ensemble simulations';
    disp(SMISMessage);
    warndlg(SMISMessage,'Warning')
    smis_ok=0; return
end

if add_diffusion==1
    SMISMessage='Diffusion not available in ensemble simulations';
    disp(SMISMessage);
    warndlg(SMISMessage,'Warning')
    smis_ok=0; return
end

if frametime>0 && (sampling_rate(1)*frametime*1e-3~=floor(sampling_rate(1)*frametime*1e-3) || ... 
            sampling_rate(1)*frametime*1e-3 < min_sampling_points)
        
    [suggested_sampling_rate, suggested_frametime]=get_sampling_rate_suggestion(frametime, min_sampling_points);
    if isempty(suggested_frametime) % A solution was found
        SMISMessage1=['Choose sampling rate for frametime so that frametime[s]*fluorescence_sampling_rate is an integer number and there is at least ', num2str(min_sampling_points),' samples per frametime'];
        SMISMessage2=['Suggested minimum sampling rate during frametime [Hz]: ',num2str(suggested_sampling_rate)];
    else
        SMISMessage1='No suitable fluorescence sampling rate could be suggested, change frametime !';
        SMISMessage2=['Suggested frametime [ms]: ',num2str(suggested_frametime)];
    end
    SMISMessage=[SMISMessage1 newline SMISMessage2];
    disp(SMISMessage);
    warndlg(SMISMessage,'Warning')
    smis_ok=0; return    
end

if addtime>0 && (sampling_rate(end)*addtime*1e-3~=floor(sampling_rate(end)*addtime*1e-3) || ...
        sampling_rate(end)*addtime*1e-3 < min_sampling_points)
    
    [suggested_sampling_rate, suggested_addtime]=get_sampling_rate_suggestion(addtime, min_sampling_points);  
    if isempty(suggested_addtime) % A solution was found
        SMISMessage1=['Choose sampling rate for addtime so that addtime[s]*fluorescence_sampling_rate is an integer number and there is at least ', num2str(min_sampling_points),' samples per frametime'];
        SMISMessage2=['Suggested minimum sampling rate during addtime [Hz]: ',num2str(suggested_sampling_rate)];
    else
        SMISMessage1='No suitable sampling rate for addtime could be suggested, change addtime !';
        SMISMessage2=['Suggested addtime [ms]: ',num2str(suggested_addtime)];
    end
    SMISMessage=[SMISMessage1 newline SMISMessage2];
    disp(SMISMessage);
    warndlg(SMISMessage,'Warning')
    smis_ok=0; return    
end


if size(lasers,2)>n_lasers % In case too many lasers defined
    lasers(n_lasers+1:end)=[];
end

if min([lasers.power])==0 % In case some lasers have zero power
    w_lasers_ok=find([lasers.power]>0);
    if isempty(w_lasers_ok)
        SMISMessage='All lasers have zero powers ! Aborting !';
        disp(SMISMessage);
        warndlg(SMISMessage,'Warning')
        smis_ok=0; return
    end
    n_lasers=numel(w_lasers_ok);
    lasers=lasers(w_lasers_ok);
end

for i=1:n_lasers
    if numel(lasers(i).sequence)==1; lasers(i).sequence=lasers(i).sequence*ones(1,n_images); end
end

for i=1:n_lasers
    lasers(i).name=strrep(lasers(i).name,' ','_');
    if numel(lasers(i).sequence)~=n_images
        SMISMessage=['Laser sequence should correspond to the number of collected images ! Check: ''', lasers(i).name, ''' Aborting !'];
        disp(SMISMessage);
        warndlg(SMISMessage,'Warning')
        smis_ok=0; return
    end
end

%Check for TIRF and HILO
if simul_3D==0 && any([lasers.tirf]==1)
    SMISMessage1='TIRF mode only valid in 3D simulations ! Check lasers !';
    w=find([lasers.tirf]==1);
    SMISMessage2=newline;
    for i=1:numel(w)
        SMISMessage2=[SMISMessage2 lasers(w(i)).name newline]; %#ok<*AGROW>
    end
    SMISMessage=[SMISMessage1 SMISMessage2];
    disp(SMISMessage);
    warndlg(SMISMessage,'Warning')
    smis_ok=0; return
end

% Check for fluorophore paterns distribution
for i=1:n_fluorophores
    if isempty(in_images_names{i})
        SMISMessage1=['Fluorophore: ',num2str(i)];
        SMISMessage2='A virtual sample must be loaded to run a simulation !';
        SMISMessage=[SMISMessage1 newline SMISMessage2];
        disp(SMISMessage);
        warndlg(SMISMessage,'Warning')
        smis_ok=0; return
    end
    
    if isempty(sm_par(i).n_sp_fr)
        SMISMessage1=['Fluorophore: ',num2str(i)];
        SMISMessage2='Fluorophore distribution in sub patterns has not been defined !';
        SMISMessage=[SMISMessage1 newline SMISMessage2];
        disp(SMISMessage);
        warndlg(SMISMessage,'Warning')
        smis_ok=0; return
    end
    
    if any(sm_par(i).n_mol.*sm_par(i).n_sp_fr~=floor(sm_par(i).n_mol.*sm_par(i).n_sp_fr))
        SMISMessage1=['Fluorophore: ',num2str(i)];
        SMISMessage2='Molecules cannot be properly distributed, check number of molecules or fractions in patterns !';
        SMISMessage=[SMISMessage1 newline SMISMessage2];
        disp(SMISMessage);
        warndlg(SMISMessage,'Warning')
        smis_ok=0; return
    end
end


%% Check if files already exists
for i=1:n_fluorophores
    ok=check_input(char(in_images_dir(i)),char(in_images_names(i)));
    if ok==0
        MyFile=fullfile(char(in_images_dir(i)),char(in_images_names(i)));
        SMISMessage=['Pattern file: ', MyFile, ' does not exist !'];
        disp(SMISMessage);
        warndlg(SMISMessage,'Warning')
        smis_ok=0; return
    end
end

for i=1:n_lasers
    if lasers(i).frap_mode==1
        SMISMessage='Ensemble simulation not available in FRAP mode !';
        disp(SMISMessage);
        warndlg(SMISMessage,'Warning')
        smis_ok=0; return
    end
end



% Check if  files should be overwritten
[simulation_files,continue_OK]=prepare_files(outfiledir,outfilename,lasers,force_overwrite);
if continue_OK==0
    smis_ok=3; return
end
if ~exist(outfiledir,'dir'); mkdir(outfiledir); end

% Open the diary
diary(simulation_files.diary_out); diary on;

%% Define general parameters for the laser
[lasers(:).frametime_duration]=deal(frametime);
[lasers(:).addtime_duration]=deal(addtime);

%% Prepare for monitoring of the photophysical status
sm_par=prepare_photophysical_status(sm_par,n_images);

%% Load fluorophores
for i=1:n_fluorophores
    sm_par(i).n_fluorescent_states=numel(sm_par(i).fluorescent_states);
    sm_par(i).n_bleached_states=numel(sm_par(i).bleached_states); % # of bleached states
    sm_par(i).n_photoactive_dark_states = numel(sm_par(i).photoactive_states)-sm_par(i).spectral_data.n_fluorescent_states; % number of dark states
    sm_par(i).photoactive_dark_states=setdiff(sm_par(i).photoactive_states,sm_par(i).fluorescent_states);
    sm_par(i).n_dark_states = sm_par(i).n_states-sm_par(i).n_fluorescent_states-sm_par(i).n_bleached_states; % number of dark states
    %
    if sm_par(i).n_photoactive_dark_states~=sm_par(i).spectral_data.n_photoactive_dark_states
        SMISMessage1=['The number of photoactive dark states does not match fluorophore setting: ', sm_par(i).fluorophore_name];
        SMISMessage2=['The number of defined photoactive dark states is: ', num2str(sm_par(i).n_photoactive_dark_states)];
        SMISMessage3=['The number of defined photoactive dark states in file: ', sm_par(i).fluorophore_name, '.mat is: ',num2str(sm_par(i).spectral_data.n_photoactive_dark_states)];
        SMISMessage4='See Matlab script: prepare_fluorophores_in_matlab_format.mat ';
        SMISMessage=[SMISMessage1 newline SMISMessage2 newline SMISMessage3 newline SMISMessage4];
        disp(SMISMessage);
        warndlg(SMISMessage,'Warning')
        smis_ok=0; return
    end
    
    %For fluorescent states in rapid equilibrium (via pH effects) get the fraction of each state
    if sm_par(i).pH_sensitivity==1
        for k=1:sm_par(i).n_fluorescent_states
            %Find the associated dark state
            sm_par(i).associated_dark_states(k)=find(sm_par(i).state_ids==sm_par(i).fluorescent_states(k) & ((1:sm_par(i).n_states)' ~= sm_par(i).fluorescent_states(k)));
            sm_par(i).fluorescent_fraction(k)=get_Hasselbach(sm_par(i).pH,sm_par(i).pKa(k),sm_par(i).Hill(k)); % Get fraction of fluorescent species (anionic state)
        end
    else
        sm_par(i).fluorescent_fraction=ones(1,sm_par(i).n_fluorescent_states); % Set fraction to 1
    end
    
    % Get the exctinction coeff at the various laser wavelength
    %     if sm_par(i).pH_sensitivity==0
    %         sm_par(i).spectral_data=get_exct_coeff_at_lambda_2(sm_par(i), n_lasers, lasers);
    %     elseif sm_par(i).pH_sensitivity==1
    %         sm_par(i).spectral_data=get_exct_coeff_at_lambda_3(sm_par(i), n_lasers, lasers);
    %     end
    
    % Get the exctinction coeff at the various laser wavelength
    % We need to handle that differently for ensemble level simulations
    % Here the states in rapid equilibrium need to be treated independently, and then re-equilibrated at every step
    sm_par(i).spectral_data=get_exct_coeff_at_lambda_2(sm_par(i), n_lasers, lasers);

    if sm_par(i).pH_sensitivity==1
        for k=1:sm_par(i).n_fluorescent_states
            for j=1:sm_par(i).n_subpop
                ens(i).sp(j).p(1,sm_par(i).associated_dark_states(k))=...
                    (1-sm_par(i).fluorescent_fraction(k))*ens(i).sp(j).p(1,sm_par(i).fluorescent_states(k)); %Set up initial state
                ens(i).sp(j).p(1,sm_par(i).fluorescent_states(k))=...
                    sm_par(i).fluorescent_fraction(k)*ens(i).sp(j).p(1,sm_par(i).fluorescent_states(k)); %Set up initial state
            end
        end
    end   
end

%% Define general parameters for the images
im_par=struct(...
    'n',[], ... % X size of pattern
    'm',[],... % Y size of pattern
    'nz',[], ... % Z size of pattern
    'raster',raster,...
    'binning', binning, ...
    'd1d2_dist',d1d2_dist, ...
    'd1d2_dist_sig', d1d2_dist_sig, ...
    'd2_constrained', d2_constrained, ...
    'obj', [], ... % objective parameters
    'det', [], ... % detector parameters
    'bg', [], ... % Background parameters
    'filters', [], ... % Filters parameters
    'n_images', n_images, ...
    'current_frame', 0, ...
    'frametime', frametime, ...
    'addtime', addtime, ...
    'during_frametime', [], ...
    'add_drift', add_drift, ...
    'add_diffusion', add_diffusion, ...
    'rbox_ch1', 0, ...
    'rbox_ch2', 0, ...
    'distribute_in_clusters',distribute_in_clusters,...
    'simul_3D',simul_3D,...
    'psf_mode',psf_mode,...
    'psf_astigmatism_x',psf_astigmatism_x, ...
    'psf_astigmatism_y',psf_astigmatism_y, ...
    'psf_astigmatism_ch1_only',psf_astigmatism_ch1_only, ...
    'psf_n_zslices', psf_n_zslices, ...
    'move_non_activated_molecules', move_non_activated_molecules, ...
    'use_diffuse_psf', use_diffuse_psf, ...
    'prevent_diffusion_out_FOV', prevent_diffusion_out_FOV, ...
    'diffuse_psf_radius', diffuse_psf_radius, ...
    'sample_zcenter', sample_zcenter, ...
    'depth_of_focus', depth_of_focus, ...
    'stochastic_spectrum_off',stochastic_spectrum_off,...
    'two_channel', two_channel, ...
    'two_channel_defocus', two_channel_defocus, ...
    'two_channel_deform', two_channel_deform, ...
    'single_CCD', single_CCD, ...
    'mol_density', 0, ...
    'fret_on', fret_on, ...
    'fluorophore_pairing_on', fluorophore_pairing_on, ...
    'apply_poisson_stat_for_lasers',apply_poisson_stat_for_lasers, ...
    'optimize_sampling_rate', optimize_sampling_rate, ...
    'minimum_oversampling', minimum_oversampling, ...
    'min_sampling_points', min_sampling_points, ...
    'debug', debug);

%% Define parameters for the single molecules, the laser and some general parameters
sms=prepare_sms(n_fluorophores, sm_par, im_par);

%Define general parameters for single molecules
for i=1:n_fluorophores
    trans_k=sm_par(i).trans_k;
    trans_q=sm_par(i).trans_q;
    sm_par(i).trans_indices_k=find(trans_k'>0); % Transitions indices for thermally activated photophysics
    sm_par(i).n_trans_k=numel(sm_par(i).trans_indices_k); % # of thermally activated transitions
    sm_par(i).trans_indices_q=find(trans_q'>0); % Transitions indices for light-activated photophysics
    sm_par(i).n_trans_q=numel(sm_par(i).trans_indices_q); % # of transitions
end

% Match pairs of fluorophores if fluorophore_pairing_on==1
if fluorophore_pairing_on==1 && n_fluorophores>1
    td_id=fluorophore_pairs(:,2);
    for i=1:numel(td_id)
        sm_par(td_id(i)).is_acceptor=1;
        sm_par(td_id(i)).td_id=fluorophore_pairs(i,1);
        sm_par(fluorophore_pairs(i,1)).td_id=td_id(i);
    end
    sms=prepare_matching(sms, sm_par, fluorophore_pairs); % Prepare the matching for pairs of fluorophores
end

% Define objective parameters
im_par.obj=struct(...
    'na', obj_na, ...
    'immersion_indice', obj_immersion_indice, ...
    'sample_indice', obj_immersion_sample, ...
    'critical_angle', asin(obj_immersion_sample/obj_immersion_indice), ...
    'mic_transmission', obj_mic_transmission, ...
    'eff_from_opening_angle', obj_eff_from_opening_angle, ...
    'eff_a',0, ...
    'eff_b',0, ...
    'mean_det_eff',obj_transmission_eff ...
    );

% Define filters parameters
im_par.filters=struct(...
    'ch1_filter', ch1_filter, ...
    'ch2_filter', ch2_filter, ...
    'exp_ch1_filter', exp_ch1_filter, ...
    'exp_ch2_filter', exp_ch2_filter, ...
    'use_exp_ch1_filter', use_exp_ch1_filter,...
    'use_exp_ch2_filter', use_exp_ch2_filter,...
    'add_dichroic', add_dichroic, ...
    'use_exp_dichroic', use_exp_dichroic,...
    'exp_dichroic_filter', exp_dichroic_filter, ...
    'dichroic_filter', dichroic_filter ...
     );

%% Create transmission patterns from filters
for i=1:n_fluorophores
    sm_par(i)=get_sm_filter_profiles(im_par,sm_par(i));
end

%% Read a high-resolution template image.

[a_h_all,n,m,nz,n_sp] = read_patterns(in_images_dir, in_images_names, im_par, simul_3D);
if isempty(a_h_all)
    SMISMessage='Virtual sample patterns could not be loaded, check filenames !';
    disp(SMISMessage);
    warndlg(SMISMessage,'Warning')
    smis_ok=0; return
end
%Update im_par
im_par.n=n; % Size of detected (binned) image
im_par.m=m;
im_par.nz=nz;
im_par.n_sp=n_sp;

if distribute_in_clusters==0
    for i=1:n_fluorophores
        sm_par(i).n_sp=numel(sm_par(i).n_sp_id);

        if numel(sm_par(i).n_sp_id)~=numel(sm_par(i).n_sp_fr)
            SMISMessage=['Fluorophore: ',sm_par(i).fluorophore_name,' not properly distributed on pattern !'];
            disp(SMISMessage);
            warndlg(SMISMessage,'Warning')
            smis_ok=0; return
        end
        
        if sm_par(i).fluorogenic==1
            if numel(sm_par(i).fluorogenicity)~=sm_par(i).n_sp
                SMISMessage=['Specified fluorogenicity must be an array of same size as # of subpatterns: ', num2str(sm_par(i).n_sp)];
                disp(SMISMessage);
                warndlg(SMISMessage,'Warning')
                smis_ok=0; return
            end
        end
        
        if n_sp(i)~=sm_par(i).n_sp
            SMISMessage1=['Incompatible # of subpatterns for fluorophore: ', num2str(i)];
            SMISMessage2=['# of specified subpatterns: ', num2str(sm_par(i).n_sp)];
            SMISMessage=[SMISMessage1 newline SMISMessage2];
            disp(SMISMessage);
            warndlg(SMISMessage,'Warning')
            smis_ok=0; return
        end
        if sum(sm_par(i).n_sp_fr)~=1
            SMISMessage=['Sum of fractional occupancies (n_sp_fr) for subpatterns for fluorophore: ', num2str(i), ' must be = 1 !'];
            disp(SMISMessage);
            warndlg(SMISMessage,'Warning')
            smis_ok=0; return
        end
    end
else
    for i=1:n_fluorophores
        sm_par(i).n_sp=1; % if distribute in clusters n_sp must be 1
        sm_par(i).n_sp_fr=1; % if distribute in clusters n_sp_fr must be 1
    end
end

%find out the density of labeling molecules
mol_density=zeros(1, n_fluorophores);
for i=1:n_fluorophores
    n_pix = length(find(a_h_all(i,:,:)>=1)); % number of pixels that can be labeled in high-res image (special case for distribute_in_clusters=1)
    mol_density(i) = 1e6*sm_par(i).n_mol/(n_pix*(raster/binning)^2); % density of mol/um^2
    theoretical_resol = 1000*2/sqrt(mol_density(i)); % factor of 2 for Nyquist
    disp(['Molecular density per micron^2 for fluorophore: ', num2str(i), ' :', num2str(mol_density(i))]);
    disp(['Theoretical achievable resolution based on Nyquist [nm]: ', num2str(theoretical_resol)]);
    disp(['Molecular density per pixel in high-resolution image: ', num2str(sm_par(i).n_mol/n_pix)]);
    disp(['Molecular density per pixel in recorded image: ', num2str(sm_par(i).n_mol/(n_pix/binning^2))]);
end
im_par.mol_density=mol_density;

%% Place single molecules on the high-resolution image

disp('Placing molecules on patterns ...')
if distribute_in_clusters==1
    if fluorophore_pairing_on==1
        SMISMessage='fluorophore_pair_on mode must be set to 0 when distribute_in_clusters=1!';
        disp(SMISMessage);
        warndlg(SMISMessage,'Warning')
        smis_ok=0; return
    end
    [sms,a_h_all_sm]=place_single_multimer_molecules(sms, a_h_all, sm_par, im_par, n_clusters);
    if isempty([sms(1).sm.x])
        SMISMessage='Molecules could not be placed: Are you sure that variable distribute_in_clusters should be set to 1 ?';
        disp(SMISMessage);
        warndlg(SMISMessage,'Warning')
        smis_ok=0; return
    end
    
else
    [sms,a_h_all_sm,sm_par,ok]=place_single_molecules(sms, a_h_all, sm_par, im_par);
    if ~ok
        smis_ok=0; return
    end
    if fluorophore_pairing_on==1 && n_fluorophores>1
        [sms,sm_par,a_h_all_sm,ok]=place_single_td_molecules(sms, a_h_all_sm, a_h_all, fluorophore_pairs, sm_par, im_par);
        if ~ok
            smis_ok=0; return
        end
    end
end
disp('Done ...')

a_h_all_sm(a_h_all_sm > 0)=100; % Reset a_h_all

%Account for limited maturation
[sms,sm_par,a_h_all_sm]=check_maturation_level(sms, sm_par, im_par, a_h_all_sm);
if fluorophore_pairing_on==1 && n_fluorophores>1
    sms=check_matching(sms, fluorophore_pairs);
end


% Get the XY positions for P
for i=1:n_fluorophores
    % TO BE CHANGED FOR SEVERAL sp
    ens(i).sp(1).x = mean([sms(i).sm.x]); % Average X position
    ens(i).sp(1).y = mean([sms(i).sm.y]); % Average Y position
    if simul_3D==1
        ens(i).sp(1).z = mean([sms(i).sm.z]); % Average Y position
    end
end


% Orient single molecules if anisotropy set to 1
sms=prepare_dipole_orientation(sms,n_fluorophores,sm_par);


% Merge all images in a_h_all to create a rgb image
if simul_3D==0
    im_par.true_im_rgb=prepare_im_rgb(a_h_all_sm, im_par);
elseif simul_3D==1
    im_par.true_3D_kernel=prepare_3D_kernel(a_h_all_sm);
end

%% Prepare the Camera (low-resolution images)

%Display image size
disp(['Imaged sample size is [um, um]: ',num2str(0.001*raster*[n,m])])
if simul_3D==1
    disp(['Sample thickness is [um]: ',num2str(0.001*raster*nz)])
end

%% Initialize lasers
% Get excitation beam profiles
lasers=get_beam_profile(n_lasers, lasers, im_par);

%% Get detector efficiency
if im_par.obj.eff_from_opening_angle==1
    opening_angle=asin(im_par.obj.na/im_par.obj.immersion_indice);
    disp(['Objective half opening angle: ',num2str(opening_angle*180/pi)]);
    [im_par.obj.eff_a, im_par.obj.eff_b, im_par.obj.mean_det_eff]=get_det_eff(opening_angle,im_par.obj.mic_transmission);
    disp(['Mean objective detection efficiency: ',num2str(im_par.obj.mean_det_eff)]);
end

%% Check sampling rate
for i=1:n_fluorophores
    sm_par(i)=check_sampling_rate(sms(i), lasers, sm_par(i), im_par);
end

%% RUN SIMULATION

disp('************ Now running simulation ... ***************');

for frame=1:n_images
    
    %Handle signal from GUI to stop simulations if not in batch mode
    if ~isempty(app)
        pause(min([exec_time*0.05,0.5])); % Necessary to catch stop signal from GUI for very fast running simulations
        if app.STOPCheckBox.Value==1
            smis_interrupt=1;
            disp('Stopping SMIS Simulation !');
            break
        end

        if 100*frame/n_images==floor(100*frame/n_images)
            app.SMISProgressGauge.Value=100*frame/n_images;
        end
    end
    
    disp(['Processing frame #: ', num2str(frame)]);
    for i=1:n_lasers
        disp(['Laser ',num2str(i),' (', num2str(lasers(i).name), ...
            '): Power: ', num2str(lasers(i).power*lasers(i).sequence(frame)/100), ...
            ' mW, Wavelength: ', num2str(lasers(i).wavelength),' nm']);
    end
    im_par.current_frame=frame;
    
    %%
    tic % Start timing the whole process
    
    %% Process photophysical changes during addtime and frametime
    % ****** This is a key script ********
    
    for i=1:n_fluorophores
        % Process ensemble photophysics to get state populations
        [ens(i), sm_par(i)]=process_ensemble_photophysics(ens(i), lasers, sm_par(i), im_par); % Do not update sm_par.sampling_rate to keep reasonnable value     
    end
    % ****** End of key script ********
       
    %% Get the detectable fluorescence (# of emitted photons per molecule) from the ensemble
    %Because the total state population is arbitrarily given the value of one, the calculated number of emitted photons corresponds to a number per molecule.
    
    ens=get_ensemble_fluorescence(ens, sm_par, im_par); % Do not update sm_par.sampling_rate to keep reasonnable value
        
    %% Check if data collection to be continued or stopped
    stop_data_collection=check_stop(sms, sm_par, im_par);
    if stop_data_collection==1
        break
    end
    
    %% Finish timing the whole process
    exec_time=toc;
    disp(['Elapsed time for this frame: ',num2str(exec_time),' seconds']);
    
end % End of frame #


%% Plot Data

plot_ensemble_data(ens,im_par,sm_par,smis_par.smis_title);


%% Save data (Even if simulation was interrupted)
save(simulation_files.data_out,'im_par','sm_par','ens','lasers'); %Save data
disp('Good, your ensemble SMIS Data are saved !');
if smis_interrupt==0   
    smis_ok=1; % If program reaches this point, SMIS simulation was successful
else
    smis_ok=2; % If program reaches this point, SMIS simulation was canceled by the SMIS GUI
end

%% Stop the diary
diary off


%% End of program
