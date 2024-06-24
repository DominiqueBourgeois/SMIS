function [sm_cell, sm_par] = process_photophysics_pct(sm_cell, lasers, sm_par, im_par)

%
% PURPOSE:
%   Main script to process photophysical state from single molecules
%
% INPUTS:
%   sm_cell: the single molecules (with coordinates on high-resolution image) in raster units
%   lasers: the lasers
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   sm_cell: the single molecules updated for photophysical state
%	sm_par: the sm parameters eventually updated for sampling rate
%
% MODIFICATION HISTORY:
%	D.Bourgeois, July 2019.
%	D.Bourgeois, June 2022, adapted to parallel computing
%	D.Bourgeois, September 2022, optimized for parallel computing

%First define general parameters
trans_indices_k=sm_par.trans_indices_k; % Indices for thermal transitions
trans_indices_q=sm_par.trans_indices_q; % Indices for photo activated transitions
trans_indices = union(trans_indices_k,trans_indices_q); % union between these two

n_trans=numel(trans_indices); % overall number of transitions to consider
w_events(1:n_trans)=struct('w', zeros(1,sm_par.n_states)); % initialize array of indices were transitions might be detected

starting_states=zeros(2,numel(trans_indices));
ending_states=zeros(2,numel(trans_indices));
starting_states(1,:)=trans_indices;
ending_states(1,:)=trans_indices;
[ending_states(2,:),starting_states(2,:)] = ind2sub(size(sm_par.trans_k),trans_indices); % Pour chaque transition, état de départ !

% Each mol will be processed individually (eventual FRET will be considered later)
addtime=im_par.addtime;

%Extract the needed fields for sm_cell
sm=sm_cell(sm_par.idx.process_photophysics.indices,:);

% tic
parfor (k=1:sm_par.n_mol_eff, im_par.parforArg)
%for k=1:sm_par.n_mol_eff
    if addtime>0 %Process photophysics during addtime
        during_frametime=0;
        sm(:,k)=get_state_evolution_pct(sm(:,k),lasers, sm_par, im_par, w_events, starting_states, ending_states, during_frametime);
    end
    %Now process photophysics during frametime
    during_frametime=1;
    sm(:,k)=get_state_evolution_pct(sm(:,k), lasers, sm_par, im_par, w_events, starting_states, ending_states, during_frametime);
end
% toc

%Fill up sm_cell
sm_cell(sm_par.idx.process_photophysics.indices,:)=sm;



