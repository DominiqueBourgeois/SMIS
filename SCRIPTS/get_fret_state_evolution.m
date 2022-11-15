function [donor, acceptor]=get_fret_state_evolution(donor, donor_par, acceptor, acceptor_par, ...
    lasers, im_par, w_events_D, starting_states_D, ending_states_D, ...
    w_events_A, starting_states_A, ending_states_A, state_link_D, state_link_A, fluo_states_D, fluo_states_A, ...
    R0_link_D, R0_link_A, during_frametime)

%
% PURPOSE:
%   Get the evolution of photophysical state of an acceptor molecule in
%   FRET
%
% INPUTS:
%   acceptor: the acceptor molecule (with coordinates on high-resolution image in raster units)
%	acceptor_par: the acceptor parameters
%   donor: the donor molecule (with coordinates on high-resolution image in raster units)
%	donor_par: the donor parameters
%   lasers: the lasers
%	im_par: the imaging parameters
%   w_events_D: structure array of indices were donor transitions might be detected (only the variable is input, its content is derived in the present routine)
%   starting_states_D: array where starting states for the different donor transitions are specified
%   ending_states_D: array where ending states for the different donor transitions are specified
%   w_events_A: structure array of indices were acceptor transitions might be detected (only the variable is input, its content is derived in the present routine)
%   starting_states_A: array where starting states for the different acceptor transitions are specified
%   ending_states_A: array where ending states for the different acceptor transitions are specified
%   state_link_D: structure array equivalent to starting_states_D and ending_states_D, but for submatrix of trans_k and trans_q with defined starting states
%   state_link_A: structure array equivalent to starting_states_A and ending_states_A, but for submatrix of trans_k and trans_q with defined starting states
%   fluo_states_D: the active fluorescent states of the donor
%   fluo_states_A: the active fluorescent states of the acceptor
%   R0_link_D: structure with ids of fluorescent states and associated dark states in rapid equilibrium for donor
%   R0_link_A: structure with ids of fluorescent states and associated dark states in rapid equilibrium for acceptor
%   during_frametime: set to 1 if state evolution to be processed during frametime, 0 if during addtime
%
% OUTPUTS:
%   donor: the donor molecule updated for energy transfer in each fluorescent state
%   acceptor: the acceptor molecule updated for photophysical state
%
% MODIFICATION HISTORY:
%	D.Bourgeois, July 2019.
%	D.Bourgeois, June 2020. Update state_trace: replace state_id by state
%	for better monitoring
%	D.Bourgeois, July 2020, full FRET processing
%	D.Bourgeois, June 2022, corrected bug for donor bleached or acceptor
%	bleached case was running twice through addtime + frametime (both in
%	this routine and in process_fret_photophysics !)
%	D.Bourgeois, June 2022, adapted to parallel computing


plot_traces=0; % Set to 1 if every trace to be displayed
plot_fluo_trace=0; % Set to 1 to see fluorescent state appears

%%   if the donor is bleached: process acceptor only
if donor.bleached>0
    %Reset state trace
    donor.state_trace=[];
    donor.fluo_trace=[];
    %Reset fret photons
    donor.given_fret_photons=zeros(1,donor_par.n_fluorescent_states);
    donor.received_fret_photons=zeros(1,donor_par.n_fluorescent_states);
    acceptor.given_fret_photons=zeros(1,acceptor_par.n_fluorescent_states);
    acceptor.received_fret_photons=zeros(1,acceptor_par.n_fluorescent_states);
  
    acceptor=get_state_evolution(acceptor,lasers, acceptor_par, im_par, w_events_A, starting_states_A, ending_states_A, during_frametime);
  
    return
end

%%   if the acceptor is bleached: process donor only
if acceptor.bleached>0
    %Reset state trace
    acceptor.state_trace=[];
    acceptor.fluo_trace=[];
    %Reset fret photons
    donor.given_fret_photons=zeros(1,donor_par.n_fluorescent_states);
    donor.received_fret_photons=zeros(1,donor_par.n_fluorescent_states);
    acceptor.given_fret_photons=zeros(1,acceptor_par.n_fluorescent_states);
    acceptor.received_fret_photons=zeros(1,acceptor_par.n_fluorescent_states);

    donor=get_state_evolution(donor,lasers, donor_par, im_par, w_events_D, starting_states_D, ending_states_D, during_frametime);

    return
end

%% Run general case
im_par.during_frametime=during_frametime; 
start_time=0; % Starting time for trace

%Reset FRET photons
donor.given_fret_photons=zeros(1,donor_par.n_fluorescent_states);
donor.received_fret_photons=zeros(1,donor_par.n_fluorescent_states);
acceptor.given_fret_photons=zeros(1,acceptor_par.n_fluorescent_states);
acceptor.received_fret_photons=zeros(1,acceptor_par.n_fluorescent_states);

% Set sampling rates for donor and acceptor
if im_par.during_frametime==1
    end_time=1e-3*im_par.frametime; % Ending time for trace
    
    %Depending wether sampling_rate has already been calculated choose initial sampling rate
    if ~isnan(donor.sampling_rate(1))
        sampling_rate_D=donor.sampling_rate(1);
    else
        sampling_rate_D=donor_par.start_sampling_rate(1);
    end
    if ~isnan(acceptor.sampling_rate(1))
        sampling_rate_A=acceptor.sampling_rate(1);
    else
        sampling_rate_A=acceptor_par.start_sampling_rate(1);
    end
else
    end_time=1e-3*im_par.addtime;
    
    %Depending wether sampling_rate has already been calculated choose initial sampling rate
    if ~isnan(donor.sampling_rate(2))
        sampling_rate_D=donor.sampling_rate(2);
    else
        sampling_rate_D=donor_par.start_sampling_rate(2);
    end
    if ~isnan(acceptor.sampling_rate(2))
        sampling_rate_A=acceptor.sampling_rate(2);
    else
        sampling_rate_A=acceptor_par.start_sampling_rate(2);
    end
end

% Get # of absorbed photons per second for each donor and acceptor states
Ns_D_lasers=sampling_rate_D*get_number_of_absorbed_photons(donor, lasers, donor_par, im_par, sampling_rate_D);
Ns_A_lasers=sampling_rate_A*get_number_of_absorbed_photons(acceptor, lasers, acceptor_par, im_par, sampling_rate_A);

% FRET efficiency for the pair
DA_distance=10*sqrt((donor.x-acceptor.x)^2+(donor.y-acceptor.y)^2+(donor.z-acceptor.z)^2)*im_par.raster/im_par.binning; % Distance in A
E_D=1./(1+(DA_distance./donor_par.R0_D).^6); % Fret efficiency from donor to acceptor
E_A=1./(1+(DA_distance./acceptor_par.R0_D).^6); % Fret efficiency from acceptor to donor

%Eventually correct for fraction of population in fluorescent or dark
%states in rapid equilibrum
if donor_par.pH_sensitivity==1
    E_A(:,R0_link_D.w_fluo_state)=donor_par.fluorescent_fraction(R0_link_D.id_fluo_state).*E_A(:,R0_link_D.w_fluo_state);
    E_A(:,R0_link_D.w_associated_dark)=(1-donor_par.fluorescent_fraction(R0_link_D.id_associated_dark)).*E_A(:,R0_link_D.w_associated_dark);
end

%Eventually correct for fraction of population in fluorescent or dark
%states in rapid equilibrum
if acceptor_par.pH_sensitivity==1
    E_D(:,R0_link_A.w_fluo_state)=acceptor_par.fluorescent_fraction(R0_link_A.id_fluo_state).*E_D(:,R0_link_A.w_fluo_state);
    E_D(:,R0_link_A.w_associated_dark)=(1-acceptor_par.fluorescent_fraction(R0_link_A.id_associated_dark)).*E_D(:,R0_link_A.w_associated_dark);
end


%Reset state trace
donor.state_trace=[];
acceptor.state_trace=[];

% Loop to construct the traces
while start_time<end_time
    % Is the donor in a fluorescent (FRET) state and acceptor is
    % photosensitive
    %     if any(donor.state==donor_par.fluorescent_states) && any(acceptor.state==donor_par.state_ids(donor_par.R0_D_index)) % Careful, we have to take species in rapid equilirium into account
    % Corrected 23/03/2021
    if any(donor.state==donor_par.fluorescent_states) && any(acceptor.state==acceptor_par.state_ids(acceptor_par.R0_A_index)) % Careful, we have to take species in rapid equilirium into account
        %         N_FRET_A_received=E_D(donor.state==fluo_states_D,acceptor.state==donor_par.state_ids(donor_par.R0_D_index)).*Ns_D_lasers(donor.state); %Compute the # of photons received by acceptor per second
        % Corrected 23/03/2021
        N_FRET_A_received=E_D(donor.state==fluo_states_D,acceptor.state==acceptor_par.state_ids(acceptor_par.R0_A_index)).*Ns_D_lasers(donor.state); %Compute the # of photons received by acceptor per second
        N_FRET_D_given=sum(N_FRET_A_received); % Total # of photons given by donor in fluorescent state per second
        % N_FRET_D_given is the # of photon given by fluorescent states only,
        % but its dimension should be 2 is fluorescent state is in equilibrium with a dark state
        if donor_par.pH_sensitivity==1
            id_dark=donor_par.associated_dark_states(donor.state==donor_par.fluorescent_states);
            id_fluo=donor.state;
            if id_fluo<id_dark
                id_D=1;
            else
                id_D=2;
            end
            N_FRET_D_given=N_FRET_D_given*[id_fluo<id_dark, id_fluo>id_dark];
        end
    else % No FRET
        N_FRET_A_received=0;
        N_FRET_D_given=0; % No FRET photon given by donor
        id_D=1;
    end
    
    % Is the acceptor in a fluorescent (FRET) state and donor is
    % photosensitive
    %     if any(acceptor.state==acceptor_par.fluorescent_states) && any(donor.state==acceptor_par.state_ids(acceptor_par.R0_D_index)) % Careful, we have to take species in rapid equilirium into account
    % Error corrected 23/03/2021
    if any(acceptor.state==acceptor_par.fluorescent_states) && any(donor.state==donor_par.state_ids(donor_par.R0_A_index)) % Careful, we have to take species in rapid equilirium into account
        %         N_FRET_D_received=E_A(acceptor.state==fluo_states_A,donor.state==acceptor_par.state_ids(acceptor_par.R0_D_index))*Ns_A_lasers(acceptor.state); %Compute the # of photons received by donor per second
        % Error corrected 23/03/2021
        N_FRET_D_received=E_A(acceptor.state==fluo_states_A,donor.state==donor_par.state_ids(donor_par.R0_A_index))*Ns_A_lasers(acceptor.state); %Compute the # of photons received by donor per second
        N_FRET_A_given=sum(N_FRET_D_received); % Total # of photons given by acceptor in fluorescent state per second
        % N_FRET_A_given is the # of photon given by fluorescent states only,
        % but its dimension should be 2 is fluorescent state is in equilibrium with a dark state
        if acceptor_par.pH_sensitivity==1
            id_dark=acceptor_par.associated_dark_states(acceptor.state==acceptor_par.fluorescent_states);
            id_fluo=acceptor.state;
            if id_fluo<id_dark
                id_A=1;
            else
                id_A=2;
            end
            N_FRET_A_given=N_FRET_A_given*[id_fluo<id_dark, id_fluo>id_dark];
        end
    else % No FRET
        N_FRET_D_received=0;
        N_FRET_A_given=0; % No FRET photon given by donor
        id_A=1;
    end
    
    % Get the # of photons absorbed by D and A per sampling time
    N_D=(Ns_D_lasers(donor_par.state_ids==donor.state)-N_FRET_D_given+N_FRET_D_received)/sampling_rate_D;
    N_A=(Ns_A_lasers(acceptor_par.state_ids==acceptor.state)-N_FRET_A_given+N_FRET_A_received)/sampling_rate_A;
    
    % Compute the subtraces until an eventual first event happens
    [subtrace_D,sampling_rate_D]=get_fret_subtrace(donor,donor_par,N_D,sampling_rate_D, start_time,end_time, im_par,w_events_D, state_link_D);
    [subtrace_A,sampling_rate_A]=get_fret_subtrace(acceptor,acceptor_par,N_A,sampling_rate_A, start_time,end_time, im_par,w_events_A, state_link_A);
    
    %Look at the time of the first event for both D and A and update the subtraces
    
    if subtrace_D(1,end)<subtrace_A(1,end) % There was an event first for D, so update A and D state
        subtrace_A=[subtrace_A(:,1),[subtrace_D(1,end);subtrace_A(2,1)]];
        start_time=subtrace_D(1,end); % New start time
        donor.state=subtrace_D(2,end);
        donor.state_id = donor_par.state_ids(donor.state); % update sm state_id
    elseif subtrace_A(1,end)<subtrace_D(1,end) % There was an event first for A, so update D and A state
        subtrace_D=[subtrace_D(:,1),[subtrace_A(1,end);subtrace_D(2,1)]];
        start_time=subtrace_A(1,end); % New start time
        acceptor.state=subtrace_A(2,end);
        acceptor.state_id = acceptor_par.state_ids(acceptor.state); % update sm state_id
    else % Nothing has occurred or possibly a state transition at exactly the same time
        start_time=subtrace_A(1,end); % New start time
        donor.state=subtrace_D(2,end);
        acceptor.state=subtrace_A(2,end);
        donor.state_id = donor_par.state_ids(donor.state); % update sm state_id
        acceptor.state_id = acceptor_par.state_ids(acceptor.state); % update sm state_id
    end
    
    %If donor started subtrace in a fluorescence state, get the number of
    %photons received (as acceptor) or given (as donor)
    if any(subtrace_D(2,1),donor_par.fluorescent_states)
        % Length of subtrace
        subtrace_T=start_time-subtrace_D(1,1); % end of subtrace is now = start_time and beginning = subtrace_D(1,1)
        if donor_par.pH_sensitivity==1
            %id_D is the index (defined above of fluorescent state in
            %N_FRET_D_given and N_FRET_D_received
            donor.given_fret_photons(subtrace_D(2,1)==donor_par.fluorescent_states)=donor.given_fret_photons(subtrace_D(2,1)==donor_par.fluorescent_states)+subtrace_T*N_FRET_D_given(id_D);
            donor.received_fret_photons(subtrace_D(2,1)==donor_par.fluorescent_states)=donor.received_fret_photons(subtrace_D(2,1)==donor_par.fluorescent_states)+subtrace_T*N_FRET_D_received(id_D);
        else
            donor.given_fret_photons(subtrace_D(2,1)==donor_par.fluorescent_states)=donor.given_fret_photons(subtrace_D(2,1)==donor_par.fluorescent_states)+subtrace_T*N_FRET_D_given;
            donor.received_fret_photons(subtrace_D(2,1)==donor_par.fluorescent_states)=donor.received_fret_photons(subtrace_D(2,1)==donor_par.fluorescent_states)+subtrace_T*N_FRET_D_received;
        end
    end
    
    %If acceptor started subtrace in a fluorescence state, get the number of
    %photons received (as acceptor) or given (as donor)
    if any(subtrace_A(2,1),acceptor_par.fluorescent_states)
        % Length of subtrace
        subtrace_T=start_time-subtrace_A(1,1); % end of subtrace is now = start_time and beginning = subtrace_A(1,1)
        if acceptor_par.pH_sensitivity==1
            %id_A is the index (defined above of fluorescent state in
            %N_FRET_A_given and N_FRET_A_received
            acceptor.given_fret_photons(subtrace_A(2,1)==acceptor_par.fluorescent_states)=acceptor.given_fret_photons(subtrace_A(2,1)==acceptor_par.fluorescent_states)+subtrace_T*N_FRET_A_given(id_A);
            acceptor.received_fret_photons(subtrace_A(2,1)==acceptor_par.fluorescent_states)=acceptor.received_fret_photons(subtrace_A(2,1)==acceptor_par.fluorescent_states)+subtrace_T*N_FRET_A_received(id_A);
        else
            acceptor.given_fret_photons(subtrace_A(2,1)==acceptor_par.fluorescent_states)=acceptor.given_fret_photons(subtrace_A(2,1)==acceptor_par.fluorescent_states)+subtrace_T*N_FRET_A_given;
            acceptor.received_fret_photons(subtrace_A(2,1)==acceptor_par.fluorescent_states)=acceptor.received_fret_photons(subtrace_A(2,1)==acceptor_par.fluorescent_states)+subtrace_T*N_FRET_A_received;
        end
    end
    
    %Handle photoconversion
    if (donor.activated==0 && donor.state>=donor_par.converted_state); donor.activated=im_par.current_frame; end % update activation state if photoconversion occurred
    if (acceptor.activated==0 && acceptor.state>=acceptor_par.converted_state); acceptor.activated=im_par.current_frame; end % update activation state if photoconversion occurred
    
    %Merge the D and A subtraces
    donor.state_trace=merge_subtrace(donor.state_trace, subtrace_D);
    acceptor.state_trace=merge_subtrace(acceptor.state_trace, subtrace_A);
    
    %Handle possible bleaching
    if min(abs(donor.state-[donor_par.bleached_states]))==0
        donor.bleached=im_par.current_frame;
        donor.blinked=0;
    end
    if min(abs(acceptor.state-[acceptor_par.bleached_states]))==0
        acceptor.bleached=im_par.current_frame;
        acceptor.blinked=0;
    end
    
    %Update sampling rates
    if im_par.during_frametime==1
        donor.sampling_rate(end)=sampling_rate_D;
        acceptor.sampling_rate(end)=sampling_rate_A;
    else
        donor.sampling_rate(1)=sampling_rate_D;
        acceptor.sampling_rate(1)=sampling_rate_A;
    end
end

% define fluorescent traces
fluo_trace_D=donor.state_trace;
fluo_trace_D(2,:)=donor_par.state_ids(fluo_trace_D(2,:)); % Replace by state_ids for rapidly exchanging states
fluo_trace_D(2,~ismember(fluo_trace_D(2,:),[donor_par.fluorescent_states]))=0;
donor.fluo_trace=fluo_trace_D;  % updated fluorescent state

fluo_trace_A=acceptor.state_trace;
fluo_trace_A(2,:)=acceptor_par.state_ids(fluo_trace_A(2,:)); % Replace by state_ids for rapidly exchanging states
fluo_trace_A(2,~ismember(fluo_trace_A(2,:),[acceptor_par.fluorescent_states]))=0;
acceptor.fluo_trace=fluo_trace_A;  % updated fluorescent state

%update blinking status at the end of the traces
if min(abs(donor_par.state_ids(donor.state_trace(2,end))-[donor_par.fluorescent_states]))~=0 && donor.bleached==0
    donor.blinked=1;
else
    donor.blinked=0;
end

if min(abs(acceptor_par.state_ids(acceptor.state_trace(2,end))-[acceptor_par.fluorescent_states]))~=0 && acceptor.bleached==0
    acceptor.blinked=1;
else
    acceptor.blinked=0;
end

%update sampling rate in _par and Ns
if im_par.during_frametime==1
    donor.Ns(end,:)=Ns_D_lasers/sampling_rate_D; % Ns_D_lasers is in # photons per second; Ns is per sampling time
    acceptor.Ns(end,:)=Ns_A_lasers/sampling_rate_A;
    %donor.sampling_rate(1)=sampling_rate_D; % Already set in while loop 
    %acceptor.sampling_rate(1)=sampling_rate_A; % Already set in while loop 
else
    donor.Ns(1,:)=Ns_D_lasers/sampling_rate_D; % Ns_D_lasers is in # photons per second; Ns is per sampling time
    acceptor.Ns(1,:)=Ns_A_lasers/sampling_rate_A;
    %donor.sampling_rate(2)=sampling_rate_D; % Already set in while loop 
    %acceptor.sampling_rate(2)=sampling_rate_A; % Already set in while loop 
end

%update total traces
donor.tot_state_trace=update_tot_trace(donor, im_par);
acceptor.tot_state_trace=update_tot_trace(acceptor, im_par);

if plot_traces==1
    figure(1); clf
    subplot(2,1,1)
    if plot_fluo_trace==0
        p=plot(donor.state_trace(1,:),donor.state_trace(2,:));
    else
        p=plot(donor.fluo_trace(1,:),donor.fluo_trace(2,:));
    end
    set(p,'LineWidth',1.5);
    set(p,'Color','blue');
    ylim([0 max(donor.state_trace(2,:))]);
    xlabel('Time [s]','fontsize',8,'fontweight','b')
    ylabel('State [AU]','fontsize',8,'fontweight','b')
    title(gca,'Donor Single Molecule Trace','FontWeight','bold');
    disp(['Donor id:' , num2str(donor.id)]);
    disp(['Photophysical state at the end of trace:' , num2str(donor.state)]);
    subplot(2,1,2)
    if plot_fluo_trace==0
        p=plot(acceptor.state_trace(1,:),acceptor.state_trace(2,:));
    else
        p=plot(acceptor.fluo_trace(1,:),acceptor.fluo_trace(2,:));
    end
    set(p,'LineWidth',1.5);
    set(p,'Color','blue');
    ylim([0 max(acceptor.state_trace(2,:))]);
    xlabel('Time [s]','fontsize',8,'fontweight','b')
    ylabel('State [AU]','fontsize',8,'fontweight','b')
    title(gca,'Acceptor Single Molecule Trace','FontWeight','bold');
    disp(['Acceptor id:' , num2str(acceptor.id)]);
    disp(['Photophysical state at end of trace:' , num2str(acceptor.state)]);
    drawnow
    input('Ok ? ','s');
end






