function [donor_cell, donor_par, acceptor_cell, acceptor_par] = process_fret_photophysics_pct(donor_cell, donor_par, acceptor_cell, acceptor_par, lasers, im_par)

%
% PURPOSE:
%   Main script to process photophysical state from FRET pair. Note that
%   FRET will be processed in both ways: from donor to acceptor, but also
%   from acceptor to donor.
%
% INPUTS:
%   donor_cell: the donor single molecules 
%	donor_par: the sm parameters of the donor
%   acceptor_cell: the acceptor single molecules
%	acceptor_par: the sm parameters of the acceptor
%   lasers: the lasers
%	im_par: the imaging parameters
%
% OUTPUTS:
%   donor: the donor single molecules updated for photophysical state
%	donor_par: the updated sm parameters of the donor
%   acceptor: the acceptor single molecules updated for photophysical state
%	acceptor_par: the updated sm parameters of the acceptor
%
% MODIFICATION HISTORY:
%	D.Bourgeois, July 2020.
%	D.Bourgeois, September 2022, optimized for parallel computing


%% Extract the needed fields for donor and acceptor
donor=donor_cell(donor_par.idx.process_fret_photophysics.indices,:);
acceptor=acceptor_cell(donor_par.idx.process_fret_photophysics.indices,:);

% Get the proper indices for donor and acceptor
id_idx=4;
matched_idx=19;

%% Define general parameters for donor
%This is the general starting_states and ending_states when all transitions
%are treated at once, which will only be the case for non paired dyes
trans_indices_k_D=donor_par.trans_indices_k; % Indices for thermal transitions
trans_indices_q_D=donor_par.trans_indices_q; % Indices for photo activated transitions
trans_indices_D = union(trans_indices_k_D,trans_indices_q_D); % union between these two

n_trans_D=numel(trans_indices_D); % overall number of transitions to consider
w_events_D(1:n_trans_D)=struct('w', zeros(1,donor_par.n_states)); % initialize array of indices were transitions might be detected

starting_states_D=zeros(2,numel(trans_indices_D));
ending_states_D=zeros(2,numel(trans_indices_D));
starting_states_D(1,:)=trans_indices_D;
ending_states_D(1,:)=trans_indices_D;
[ending_states_D(2,:),starting_states_D(2,:)] = ind2sub(size(donor_par.trans_k),trans_indices_D); % Pour chaque transition, état de départ !

% for paired dyes we define starting and ending states differently to adapt to the fact that only one transition will be searched at each subtrace
state_link_D(1:donor_par.n_states)=struct('index',[],'starting_states',[],'ending_states',[]);
for i=1:donor_par.n_states % look which transition can occur for each starting state
    state_link_D(i).index=find(donor_par.trans_k(donor_par.state_ids==i,:)' > 0 | donor_par.trans_q(donor_par.state_ids==i,:)'>0);
    state_link_D(i).starting_states=starting_states_D(2,donor_par.state_ids(starting_states_D(2,:))==i);
    state_link_D(i).ending_states=ending_states_D(2,donor_par.state_ids(starting_states_D(2,:))==i);
end

%Define in R0_index which are fluorescent states and associated dark states
R0_link_D=struct('w_fluo_state',[],'id_fluo_state',[],'w_associated_dark',[],'id_associated_dark',[]);
if donor_par.pH_sensitivity==1
    %Which acceptor state is fluorescent ?
    [tmp_1,tmp_2]=ismember(acceptor_par.R0_D_index,donor_par.fluorescent_states);
    R0_link_D.w_fluo_state=find(tmp_1==1); % Index of fluorescent states in R0_D_index
    R0_link_D.id_fluo_state=tmp_2(R0_link_D.w_fluo_state); % Index of  fluorescent state in donor_par.fluorescent_states
    %find the associated dark states
    [tmp_1,tmp_2]=ismember(acceptor_par.R0_D_index,donor_par.associated_dark_states(R0_link_D.id_fluo_state));
    R0_link_D.w_associated_dark=find(tmp_1==1);
    R0_link_D.id_associated_dark=tmp_2(R0_link_D.w_associated_dark); % Which fluorescent state id
end


%% Define general parameters for acceptor
trans_indices_k_A=acceptor_par.trans_indices_k; % Indices for thermal transitions
trans_indices_q_A=acceptor_par.trans_indices_q; % Indices for photo activated transitions
trans_indices_A = union(trans_indices_k_A,trans_indices_q_A); % union between these two

n_trans_A=numel(trans_indices_A); % overall number of transitions to consider
w_events_A(1:n_trans_A)=struct('w', zeros(1,acceptor_par.n_states)); % initialize array of indices were transitions might be detected

starting_states_A=zeros(2,numel(trans_indices_A));
ending_states_A=zeros(2,numel(trans_indices_A));
starting_states_A(1,:)=trans_indices_A;
ending_states_A(1,:)=trans_indices_A;
[ending_states_A(2,:),starting_states_A(2,:)] = ind2sub(size(acceptor_par.trans_k),trans_indices_A); % Pour chaque transition, état de départ !

% for paired dyes we define starting and ending states differently to adapt to the fact that only one transition will be searched at each subtrace
state_link_A(1:acceptor_par.n_states)=struct('index',[],'starting_states',[],'ending_states',[]);
for i=1:acceptor_par.n_states % look which transition can occur for each starting state
    state_link_A(i).index=find(acceptor_par.trans_k(acceptor_par.state_ids==i,:)' > 0 | acceptor_par.trans_q(acceptor_par.state_ids==i,:)'>0);
    state_link_A(i).starting_states=starting_states_A(2,acceptor_par.state_ids(starting_states_A(2,:))==i);
    state_link_A(i).ending_states=ending_states_A(2,acceptor_par.state_ids(starting_states_A(2,:))==i);
end

R0_link_A=struct('w_fluo_state',[],'id_fluo_state',[],'w_associated_dark',[],'id_associated_dark',[]);
if acceptor_par.pH_sensitivity==1
    %Which acceptor state is fluorescent ?
    [tmp_1,tmp_2]=ismember(donor_par.R0_D_index,acceptor_par.fluorescent_states);
    R0_link_A.w_fluo_state=find(tmp_1==1); % Index of fluorescent states in R0_D_index
    R0_link_A.id_fluo_state=tmp_2(R0_link_A.w_fluo_state); % Index of  fluorescent state in donor_par.fluorescent_states
    %find the associated dark states
    [tmp_1,tmp_2]=ismember(donor_par.R0_D_index,acceptor_par.associated_dark_states(R0_link_A.id_fluo_state));
    R0_link_A.w_associated_dark=find(tmp_1==1);
    R0_link_A.id_associated_dark=tmp_2(R0_link_A.w_associated_dark); % Which fluorescent state id
end


%% Process the eventual left over (unpaired) molecules and get indices of paired molecules
n_pairs=min([acceptor_par.n_mol_eff,donor_par.n_mol_eff]);

if acceptor_par.n_mol_eff>donor_par.n_mol_eff
    %identify the left over acceptors
    unpaired_indices_A=find([acceptor{matched_idx,:}]==0); % Indices of acceptors that are unpaired
    unpaired_indices_D=[]; % All donors are effectively paired
    paired_indices_D=1:donor_par.n_mol_eff; % All donors are effectively paired
elseif donor_par.n_mol_eff>acceptor_par.n_mol_eff
    %identify the left over donors
    unpaired_indices_A=[]; % All acceptors are effectively paired
    unpaired_indices_D=find([donor{matched_idx,:}]==0); % Indices of donors that are unpaired
    paired_indices_D=find([donor{matched_idx,:}]>0); % Indices of donors that are paired
else
    unpaired_indices_A=[]; % All acceptors are effectively paired
    unpaired_indices_D=[]; % All donors are effectively paired
    paired_indices_D=1:donor_par.n_mol_eff;
end

% Get paired_indices_A (indices of acceptors that are effectively paired
paired_indices_A=nan(1,n_pairs);
acceptor_ids=[acceptor{id_idx,:}];
for k=1:n_pairs
    paired_indice=find(acceptor_ids==donor{matched_idx,paired_indices_D(k)},1);
    if ~isempty(paired_indice)
        paired_indices_A(k)=paired_indice;
    end
end
paired_indices_A=paired_indices_A(~isnan(paired_indices_A));


%% Get fluorescent states of both donor and acceptor
fluo_states_D=donor_par.fluorescent_states(donor_par.fluorescent_states>=donor_par.initial_fluo_state);
fluo_states_A=acceptor_par.fluorescent_states(acceptor_par.fluorescent_states>=acceptor_par.initial_fluo_state);


%%  Extract the paired or impaired donor and acceptor molecules

donor_sm=donor(:,paired_indices_D);
if ~isempty(unpaired_indices_D)
    unpaired_donor_sm=donor(:,unpaired_indices_D);
end

acceptor_sm=acceptor(:,paired_indices_A);
if ~isempty(unpaired_indices_A)
    unpaired_acceptor_sm=acceptor(:,unpaired_indices_A);
end

%% Process photophysics of freting pairs

addtime=im_par.addtime; % Needed by the parfor loop
% parfor (k=1:n_pairs, im_par.parforArg)
for k=1:n_pairs
    % for k=1:n_pairs
    if addtime>0
        %Process photophysics during addtime
        during_frametime=0;
        [donor_sm(:,k), acceptor_sm(:,k)]=...
            get_fret_state_evolution_pct(donor_sm(:,k),donor_par, acceptor_sm(:,k), acceptor_par, lasers, im_par, ...
            w_events_D, starting_states_D, ending_states_D, w_events_A, starting_states_A, ending_states_A, ...
            state_link_D, state_link_A, fluo_states_D, fluo_states_A, R0_link_D, R0_link_A,during_frametime);

    end
    %Now process photophysics during frametime
    during_frametime=1;
    [donor_sm(:,k),acceptor_sm(:,k)]=...
        get_fret_state_evolution_pct(donor_sm(:,k),donor_par, acceptor_sm(:,k), acceptor_par, lasers, im_par, ...
        w_events_D, starting_states_D, ending_states_D, w_events_A, starting_states_A, ending_states_A, ...
        state_link_D, state_link_A, fluo_states_D, fluo_states_A, R0_link_D, R0_link_A,during_frametime);

end


%% process left over donors and acceptors
% Donors
if ~isempty(unpaired_indices_D)
    parfor (k=1:numel(unpaired_indices_D), im_par.parforArg)
        %for k=1:numel(unpaired_indices_D)
        if addtime>0
            during_frametime=0;
            unpaired_donor_sm(:,k)=get_state_evolution_pct(unpaired_donor_sm(:,k),lasers, donor_par, im_par, w_events_D, starting_states_D, ending_states_D, during_frametime);
        end
        during_frametime=1;
        unpaired_donor_sm(:,k)=get_state_evolution_pct(unpaired_donor_sm(:,k),lasers, donor_par, im_par, w_events_D, starting_states_D, ending_states_D,during_frametime);
    end
end

% Acceptors
if ~isempty(unpaired_indices_A)
    parfor (k=1:numel(unpaired_indices_A), im_par.parforArg)
    %for k=1:numel(unpaired_indices_A)
        if addtime>0
            during_frametime=0;
            unpaired_acceptor_sm(:,k)=get_state_evolution_pct(unpaired_acceptor_sm(:,k),lasers, acceptor_par, im_par, w_events_A, starting_states_A, ending_states_A, during_frametime);
        end
        during_frametime=1;
        unpaired_acceptor_sm(:,k)=get_state_evolution_pct(unpaired_acceptor_sm(:,k),lasers, acceptor_par, im_par, w_events_A, starting_states_A, ending_states_A, during_frametime);
    end
end

%% Reassign to donor.sm and acceptor.sm


donor(:,paired_indices_D)=donor_sm;
if ~isempty(unpaired_indices_D)
    donor(:,unpaired_indices_D)=unpaired_donor_sm;
end

acceptor(:,paired_indices_A)=acceptor_sm;
if ~isempty(unpaired_indices_A)
    acceptor(:,unpaired_indices_A)=unpaired_acceptor_sm;
end


%% Finally set the processing status
donor_par.processing_done=im_par.current_frame;
acceptor_par.processing_done=im_par.current_frame;

%Fill up donor_cell and acceptor_cell
donor_cell(donor_par.idx.process_fret_photophysics.indices,:)=donor;
acceptor_cell(donor_par.idx.process_fret_photophysics.indices,:)=acceptor;
