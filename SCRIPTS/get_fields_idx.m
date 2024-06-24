function idx=get_fields_idx(fn, im_par)

%
% PURPOSE:
%   Get the different field indices to enable easy working sms.sm_cell in
%   SMIS
%
% INPUTS:
%   fn: field names of sms.sm structure
%	im_par: the imaging parameters
%
% OUTPUTS:
%   idx: a structure containing the field names positions 
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2022, adapted to parallel computing

%% Extract fields for process_photophysics.m
% Get relevant indices in fn
w = matches(fn,{'x','y','z','id','theta','phi','state',...
    'state_id','state_trace','fluo_trace','tot_state_trace','Ns','sampling_rate',...
    'bleached','activated','blinked'});

idx.process_photophysics.indices=find(w==1);

%% Extract fields for process_fret_photophysics.m
if im_par.fret_on==1
    % Get relevant indices in fn
    w = matches(fn,{'x','y','z','id','theta','phi','state',...
        'state_id','state_trace','fluo_trace','tot_state_trace','Ns','sampling_rate',...
        'bleached','activated','blinked','given_fret_photons','received_fret_photons','matched'});

    idx.process_fret_photophysics.indices=find(w==1);
end

%% Extract fields for process_fluorescence.m
% Get relevant indices in fn
if im_par.fret_on==1
    w = matches(fn,{'sub_x','fluo_trace','Ns','sampling_rate', 'n_abs','n_em','tot_n_em',...
        'bleached','t_on','fr_t_on','given_fret_photons','received_fret_photons','matched'});
else
    w = matches(fn,{'sub_x','c_sp','fluo_trace','Ns','sampling_rate', 'n_abs','n_em','tot_n_em',...
    'bleached','t_on','fr_t_on'});
end

idx.process_fluorescence.indices=find(w==1);

%% Extract fields for process_emission_spectra.m
% % Get relevant indices in fn
w = matches(fn,{'n_em','n_phot_ch1','n_phot_ch2','tot_n_phot_ch1',...
    'tot_n_phot_ch2', 'bleached','em_spectrum'});

idx.process_emission_spectra.indices=find(w==1);

%% Extract fields for get_frame_image.m
% % Get relevant indices in fn
w = matches(fn,{'x','y','z','sub_x','sub_y','sub_z','c_sp','theta',...
    'n_phot_ch1','n_phot_ch2','n_phot_det_ch1','n_phot_det_ch2',...
    'tot_phot_det_ch1','tot_phot_det_ch2','frames_on_ch1','frames_on_ch2',...
    'bleached','fr_t_on','lx','ly','lz','le_set'});

idx.get_frame_image.indices=find(w==1);







