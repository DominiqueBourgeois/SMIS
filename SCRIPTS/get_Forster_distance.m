function [sm_par, ok]=get_Forster_distance(sm_par, im_par, fluorophore_pairs)
% NAME:
%	get_Forster_distance
%
% PURPOSE:
%       calculate Forster distance
% INPUTS:
%   sm_par: single molecule parameters
%
% OUTPUTS:
%	the updated single molecule parameters
%   ok = 1 if execution was okay, 0 if an error occurred
%
% MODIFICATION HISTORY:
%	D.Bourgeois, May 2019. July 2020: R0 from acceptor to donor also
%	calculated
%-

% for each dye pair we will calculate a Forster radius between donor (in each fluorescent state) and
% acceptor (in each photoactive state), and vice versa

ok=1;
n_pairs=size(fluorophore_pairs,1);

%First direction
for i=1:n_pairs
    %Define donor and acceptor
    if sm_par(fluorophore_pairs(i,1)).is_acceptor==0
        donor=sm_par(fluorophore_pairs(i,1));
        donor_id=fluorophore_pairs(i,1);
    else
        acceptor=sm_par(fluorophore_pairs(i,1));
        acceptor_id=fluorophore_pairs(i,1);
    end
    if sm_par(fluorophore_pairs(i,2)).is_acceptor==0
        donor=sm_par(fluorophore_pairs(i,2));
        donor_id=fluorophore_pairs(i,2);
    else
        acceptor=sm_par(fluorophore_pairs(i,2));
        acceptor_id=fluorophore_pairs(i,2);
    end
    
    if ~exist('donor','var') || ~exist('acceptor','var')
        disp('Donor acceptor pair does not exist !')
        errordlg('FRET RO cannot be calculated !','SMIS Error');
        uiwait
        ok=0;
        return
    end
    % The FRET donor states are the fluorescent states of the donor
    donor_states_index=find(donor.fluorescent_states>=donor.initial_fluo_state);
    donor_states_id=donor.fluorescent_states(donor_states_index); % Only consider states > initial state
    n_donor_states=numel(donor_states_id);
    
    % The FRET acceptor states are all the photoactive states of the acceptor
    acceptor_fluorescent_states_index=find(acceptor.fluorescent_states>=acceptor.initial_fluo_state);
    acceptor_fluorescent_states_id=acceptor.fluorescent_states(acceptor_fluorescent_states_index); % Only consider states > initial state
    n_acceptor_fluorescent_states=numel(acceptor_fluorescent_states_id);

    acceptor_dark_states_index=find(acceptor.photoactive_dark_states>=acceptor.initial_fluo_state);
    acceptor_dark_states_id=acceptor.photoactive_dark_states(acceptor_dark_states_index); % Only consider states > initial state
    n_acceptor_dark_states=numel(acceptor_dark_states_id);
    
    n_acceptor_states=n_acceptor_fluorescent_states+n_acceptor_dark_states;
    
    R0=zeros(n_donor_states, n_acceptor_states);
    R0_index=zeros(1, n_acceptor_states); 
    
    for j=1:n_donor_states
        em_spectrum=donor.spectral_data.em_spectra(donor_states_index(j)).s; %the emission spectrum of the donor
        QY=donor.quantum_yield(donor_states_index(j));
        for k=1:n_acceptor_fluorescent_states
            exc_spectrum=acceptor.spectral_data.exc_spectra(acceptor_fluorescent_states_index(k)).s; %the emission spectrum of the donor
            exc_eps=acceptor.spectral_data.exc_spectra(acceptor_fluorescent_states_index(k)).eps;
            %There is a trick here if fluorescent acceptor state is in
            %rapid equilibrium with a dark state: we must consider the eps
            %of the acceptor state when 100% present, and then the FRET
            %will be distributed to the states in equilibrium according to
            %their ratio. But exc_eps has been scaled by occupancy of the state, so we must correct that 
            if acceptor.pH_sensitivity==1
                exc_eps=exc_eps/acceptor.fluorescent_fraction(acceptor_fluorescent_states_index(k));
            end
                       
            R0(j, k)=get_R0(exc_spectrum, em_spectrum, exc_eps ,QY, im_par);
            if j==1; R0_index(k)=acceptor_fluorescent_states_id(k); end
            disp(['Fluorophore pair: ', num2str(donor_id),', ',num2str(acceptor_id),...
                ': R0 from fluorescent state: ',num2str(donor_states_id(j)),' to fluorescent state: ',num2str(acceptor_fluorescent_states_id(k)), ' is [A]: ',num2str(R0(j, k))]);
        end
        for k=1:n_acceptor_dark_states
            exc_spectrum=acceptor.spectral_data.dark_spectra(acceptor_dark_states_index(k)).s; %the emission spectrum of the donor
            exc_eps=acceptor.spectral_data.dark_spectra(acceptor_dark_states_index(k)).eps;
            % Same trick as above, if the dark state is in equilibrium with
            % a fluorescent state
            if acceptor.pH_sensitivity==1
                w_associated_fluo_state=find(acceptor.photoactive_dark_states(acceptor_dark_states_index(k))==acceptor.associated_dark_states);
                if ~isempty(w_associated_fluo_state)
                    exc_eps=exc_eps/(1-acceptor.fluorescent_fraction(w_associated_fluo_state));
                end
            end
            k0=k+n_acceptor_fluorescent_states;
            R0(j, k0)=get_R0(exc_spectrum, em_spectrum, exc_eps ,QY, im_par);
            if j==1; R0_index(k0)=acceptor_dark_states_id(k); end
            disp(['Fluorophore pair: ', num2str(donor_id),', ',num2str(acceptor_id),...
                ': R0 from fluorescent state: ',num2str(donor_states_id(j)),' to dark state: ',num2str(acceptor_dark_states_id(k)), ' is [A]: ',num2str(R0(j, k0))]);
        end   
    end
    sm_par(donor_id).R0_D=R0;
    sm_par(donor_id).R0_D_index=R0_index;
    sm_par(acceptor_id).R0_A_index=R0_index; % Also pass to acceptor dye
end

%Second direction
for i=1:n_pairs
    %Define donor and acceptor
    if sm_par(fluorophore_pairs(i,1)).is_acceptor==0
        acceptor=sm_par(fluorophore_pairs(i,1));
        acceptor_id=fluorophore_pairs(i,1);
    else
        donor=sm_par(fluorophore_pairs(i,1));
        donor_id=fluorophore_pairs(i,1);
    end
    if sm_par(fluorophore_pairs(i,2)).is_acceptor==0
        acceptor=sm_par(fluorophore_pairs(i,2));
        acceptor_id=fluorophore_pairs(i,2);
    else
        donor=sm_par(fluorophore_pairs(i,2));
        donor_id=fluorophore_pairs(i,2);
    end 
    
    if ~exist('donor','var') || ~exist('acceptor','var')
        disp('Acceptor Donor pair does not exist !')
        errordlg('FRET RO cannot be calculated !','SMIS Error');
        uiwait
        ok=0;
        return
    end
    
    % The FRET donor states are the fluorescent states of the donor
    donor_states_index=find(donor.fluorescent_states>=donor.initial_fluo_state);
    donor_states_id=donor.fluorescent_states(donor_states_index); % Only consider states > initial state
    n_donor_states=numel(donor_states_id);
    
    % The FRET acceptor states are all the photoactive states of the acceptor
    acceptor_fluorescent_states_index=find(acceptor.fluorescent_states>=acceptor.initial_fluo_state);
    acceptor_fluorescent_states_id=acceptor.fluorescent_states(acceptor_fluorescent_states_index); % Only consider states > initial state
    n_acceptor_fluorescent_states=numel(acceptor_fluorescent_states_id);

    acceptor_dark_states_index=find(acceptor.photoactive_dark_states>=acceptor.initial_fluo_state);
    acceptor_dark_states_id=acceptor.photoactive_dark_states(acceptor_dark_states_index); % Only consider states > initial state
    n_acceptor_dark_states=numel(acceptor_dark_states_id);
    
    n_acceptor_states=n_acceptor_fluorescent_states+n_acceptor_dark_states;
    
    R0=zeros(n_donor_states, n_acceptor_states); % Forster radii for each donor fluorescent state
    R0_index=zeros(1, n_acceptor_states);  % Indices of the acceptor states able to FRET
    
    for j=1:n_donor_states
        em_spectrum=donor.spectral_data.em_spectra(donor_states_index(j)).s; %the emission spectrum of the donor
        QY=donor.quantum_yield(donor_states_index(j));
        for k=1:n_acceptor_fluorescent_states
            exc_spectrum=acceptor.spectral_data.exc_spectra(acceptor_fluorescent_states_index(k)).s; %the emission spectrum of the donor
            exc_eps=acceptor.spectral_data.exc_spectra(acceptor_fluorescent_states_index(k)).eps;           
            %There is a trick here if fluorescent acceptor state is in
            %rapid equilibrium with a dark state: we must consider the eps
            %of the acceptor state when 100% present, and then the FRET
            %will be distributed to the states in equilibrium according to
            %their ratio. But exc_eps has been scaled by occupancy of the state, so we must correct that
            if acceptor.pH_sensitivity==1
                exc_eps=exc_eps/acceptor.fluorescent_fraction(acceptor_fluorescent_states_index(k));
            end
            
            R0(j, k)=get_R0(exc_spectrum, em_spectrum, exc_eps ,QY, im_par);
            if j==1; R0_index(k)=acceptor_fluorescent_states_id(k); end
            disp(['Fluorophore pair: ', num2str(donor_id),', ',num2str(acceptor_id),...
                ': R0 from fluorescent state: ',num2str(donor_states_id(j)),' to fluorescent state: ',num2str(acceptor_fluorescent_states_id(k)), ' is [A]: ',num2str(R0(j, k))]);
        end
        for k=1:n_acceptor_dark_states
            exc_spectrum=acceptor.spectral_data.dark_spectra(acceptor_dark_states_index(k)).s; %the emission spectrum of the donor
            exc_eps=acceptor.spectral_data.dark_spectra(acceptor_dark_states_index(k)).eps;
            % Same trick as above, if the dark state is in equilibrium with
            % a fluorescent state
            if acceptor.pH_sensitivity==1
                w_associated_fluo_state=find(acceptor.photoactive_dark_states(acceptor_dark_states_index(k))==acceptor.associated_dark_states);
                if ~isempty(w_associated_fluo_state)
                    exc_eps=exc_eps/(1-acceptor.fluorescent_fraction(w_associated_fluo_state));
                end
            end
            k0=k+n_acceptor_fluorescent_states;
            R0(j, k0)=get_R0(exc_spectrum, em_spectrum, exc_eps ,QY, im_par);
            if j==1; R0_index(k0)=acceptor_dark_states_id(k); end
            disp(['Fluorophore pair: ', num2str(donor_id),', ',num2str(acceptor_id),...
                ': R0 from fluorescent state: ',num2str(donor_states_id(j)),' to dark state: ',num2str(acceptor_dark_states_id(k)), ' is [A]: ',num2str(R0(j, k0))]);
        end   
    end
    sm_par(donor_id).R0_D=R0;
    sm_par(donor_id).R0_D_index=R0_index;
    sm_par(acceptor_id).R0_A_index=R0_index; % Also pass to acceptor dye
end

