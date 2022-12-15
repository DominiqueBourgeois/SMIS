function sm = get_state_evolution_pct(sm, lasers, sm_par, im_par, w_events, starting_states, ending_states, during_frametime)
%
% PURPOSE:
%   Get the evolution of photophysical state
%
% INPUTS:
%   sm: the single molecule (with coordinates on high-resolution image) in raster units
%   lasers: the lasers
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%   w_events: structure array of indices were transitions might be detected (only the variable is iput, its content is derived in the present routine)
%   starting_states: array where starting states for the different transitions are specified
%   ending_states: array where ending states for the different transitions are specified
%   during_frametime: set to 1 if state evolution to be processed during frametime, 0 if during addtime
%
% OUTPUTS:
%   sm: the single molecules updated for photophysical state
%
% MODIFICATION HISTORY:
%	D.Bourgeois, July 2019.
%	D.Bourgeois, June 2020. Update state_trace: replace state_id by state
%	for better monitoring
%	D.Bourgeois, December 2021. Sort allowed_trans_indices in increasing order of trans_p probabilities
%	D.Bourgeois, September 2022, adapted to parallel computing

print_message=0; % Set to 1 to see message upon sampling state change

% Get the proper indices in sm
x_idx=1;
y_idx=2;
z_idx=3;
id_idx=4;
theta_idx=5;
phi_idx=6;
state_idx=7;
state_id_idx=8;
state_trace_idx=9;
fluo_trace_idx=10;
tot_state_trace_idx=11;
Ns_idx=12;
sampling_rate_idx=13;
bleached_idx=14;
activated_idx=15;
blinked_idx=16;

%   if the molecule is bleached: nothing to do except empty the traces
if sm{bleached_idx}>0 
    sm{state_trace_idx}=[];
    sm{fluo_trace_idx}=[];
    return
end

plot_traces=0; % Set to 1 if every trace to be displayed
plot_fluo_trace=0; % Set to 1 to see fluorescent state appears

% get general parameters
im_par.during_frametime=during_frametime; 
% Time window to study in [s]
if im_par.during_frametime==1
    T=1e-3*im_par.frametime;
    %Depending wether used_sampling_rate has already been calculated choose
    %initial sampling rate
    if ~isnan(sm{sampling_rate_idx}(1))
        sampling_rate=sm{sampling_rate_idx}(1);
    else
        sampling_rate=sm_par.start_sampling_rate(1);
    end
else
    T=1e-3*im_par.addtime;
    if ~isnan(sm{sampling_rate_idx}(2))
        sampling_rate=sm{sampling_rate_idx}(2);
    else
        sampling_rate=sm_par.start_sampling_rate(2);
    end
end

minimum_oversampling=im_par.minimum_oversampling;
min_samples_per_frame=im_par.min_sampling_points;

Trace=[;]; % initialize trace

% *** Compute probability of transitions ***
% 1: For a thermal processe
% Relation between rate [s-1) and probability of transition
% k=p*n with k = rate [s-], p = probability (or yield) and n = # of trials per second = sampling rate
% so p_k = k/n = k[s-1]/S[s-1]

% 2: For a photo activated process:
% We should basically do a test for each absorbed photon
% Alternatively we speed up by testing every dt : if N is the # of
% photons absorbed in dt, and q the QY of the process, the probability that the process happens in dt is
% given by Q = q*(1 + (1-q)^1 + (1-q)^2 + ... + (1-q)^N-1) = 1 - (1-q)^N
% so it is sufficient to test the process every dt with effective
% transition probability p_q = 1 - (1-q)^N
% This is only okay if N >=1
% What if N < 1 ?
% Then 1/N is the # of microsecond "per absorbed photon"
% We split into 1 event per microsecond with reduced QY Q2
% These sub event are such that q = 1 - (1-Q2)^(1/N), so that
% Q2 = 1 - (1-q)^N , ie the formula is valid even if N < 1
% so in final for a photo activated process: p_q = 1 - (1-q)^N
% However, this formula does not work for q=1 (p_q will always be 1).
% Although this will practically never happen, it should be prevented by
% replacing q=1 by say q=1-1e-5

% 3: Finally for combined process which is either thermal or photo activated:
% the transition can go thermal with probability p_k. the probability of failing is 1-p_q
% and then there is a second option to go photoactivated with probability (1-p_k)*p_q
% (the reasoning is symmetrical)
% so in final for a combined process: p = p_q + p_k -p_q*p_k

% So we now calculate the transition probability matrix

% probability transition matrix for thermal processes, rates given in [s-1]
trans_p_k=sm_par.trans_k/sampling_rate;

% probability transition matrix for photo activated processes
% # of absorbed photons per sampling time in each starting photoactive state
Ns=get_number_of_absorbed_photons_pct(sm([x_idx, y_idx, z_idx, theta_idx, ...
    phi_idx, state_idx,Ns_idx]), lasers, sm_par, im_par, sampling_rate);

% equivalent probability transition matrix, rates given in [s-1], see discussion above
trans_p_q = 1-bsxfun(@power,1-sm_par.trans_q,Ns'); % transition matrix for photo activated processes (each row should be powered by a corresponding N)

% overall probability tansition matrix, rates given in [s-1]
trans_p=trans_p_k+trans_p_q-trans_p_k.*trans_p_q;

% to gain speed we should now define which state can be accessed from what
% state
[trans_p, allowed_n_trans, allowed_trans_indices]=get_allowed_transitions_pct(trans_p, sm{activated_idx}, sm_par);

% Check sampling rate. Sampling rate should be at least minimum_oversampling times more than actual rates.
% Actual rates [s-1] = trans_p*sampling_rate so that trans_p*sampling_rate must be < sampling rate/minimum_oversampling
% This translates into max(trans_p) < 1/minimum_oversampling. If not increase sampling_rate

if max(trans_p(:)) > 1/minimum_oversampling
    % Run the following as many times as necessary
    while max(trans_p(:)) > 1/minimum_oversampling
        old_sampling_rate=sampling_rate;
        new_sampling_rate=minimum_oversampling*old_sampling_rate*max(trans_p(:));
        
        % There must be an integer # of sampling in time window T
        required_number_of_samples=ceil(new_sampling_rate*T);
        new_sampling_rate=required_number_of_samples/T;
        sampling_rate=new_sampling_rate;
        %and redo the job
        trans_p_k=sm_par.trans_k/sampling_rate;
        Ns=Ns*old_sampling_rate/new_sampling_rate;
        trans_p_q = 1-bsxfun(@power,1-sm_par.trans_q,Ns'); % transition matrix for photo activated processes (each row should be powered by a corresponding N)
        trans_p=trans_p_k+trans_p_q-trans_p_k.*trans_p_q;
        [trans_p, allowed_n_trans, allowed_trans_indices]=get_allowed_transitions_pct(trans_p, sm{activated_idx}, sm_par);
    end
    optimized_sampling=1; % Set this flag for next optimization option below
    if print_message==1
        if im_par.during_frametime==1
            disp(['Sampling rate during frametime has to be increased to: ',num2str(sampling_rate),' Hz'])
        else
            disp(['Sampling rate during addtime has to be increased to: ',num2str(sampling_rate),' Hz'])
        end
    end
else
    optimized_sampling=0; % Set this flag for next optimization option below
end

% Eventually optimize sampling rate
if im_par.optimize_sampling_rate==1 && optimized_sampling==0
    [sampling_rate, trans_p] = optimize_sampling_rate(sampling_rate, sm_par, trans_p, minimum_oversampling, min_samples_per_frame, Ns, T);
end


% prevent transition probabilities of 1, if they happen (singular case, see
% above)
trans_p(trans_p==1)=1-1e-5;

% sort allowed_trans_indices in increasing order of trans_p probabilities
% This is important to treat the following case: at low oversampling, the most 
% likely transitions will have a probability close to one. 
% If they are treated first in the next for loop, they will almost always 
% preclude the selection of a low probability transition. 
% For example if the highest transition is 1-1e-5, and is treated first, 
% this transition will always be selected even if the possibility for a low 
% transition occurs, because at the corresponding w_events(i).w value, both tests will be positive.
[~,sorted_trans_p]=sort(trans_p(allowed_trans_indices));
allowed_trans_indices=allowed_trans_indices(sorted_trans_p);

N_trials=round(T*sampling_rate); % Number of trials: we will do a test every 1/sampling_rate [s]

Test=rand(allowed_n_trans,N_trials); % series of random numbers
for i=1:allowed_n_trans
    w_events(i).w = find(Test(i,:) <= trans_p(allowed_trans_indices(i)));
    Trace=horzcat(Trace,[w_events(i).w;allowed_trans_indices(i)*ones(1,length([w_events(i).w]))]);  % concatenate all these events
    % for each event of the trace, associate index equal to Trans_indices(i)
end

sm_trace=[0;sm{state_idx}]; % define (local) trace at start

% reorganize  trace in chronological order
if ~isempty(Trace)
    [~,ix]=sort(Trace,2);
    ix1=ix(1,:);
    Trace(:,:)=Trace(:,ix1(:));    % trace to process
    
    %process all events as a function of sm state
    for i=1:size(Trace,2)
        t=Trace(2,i); % Transition t occurs
        S=starting_states(2,starting_states(1,:)==t); % What is the starting state for this transition
        % Check if the transition can take place: state_id (and not state)
        % should be compared because 2 different photophysical states that have the
        % same state_id should evolve together. This was designed to
        % handle photoconversion which applies to the protonated state of a
        % PCFP, not to the anionic state.
        if sm{state_id_idx} == sm_par.state_ids(S)
            %Check also that a state change is unique at a given time
            if Trace(1,i)*T/N_trials ~= sm_trace(1,end)
                sm_trace = horzcat(sm_trace,[Trace(1,i)*T/N_trials;S]); % initial state
                sm{state_idx} = ending_states(2,ending_states(1,:)==t); % update sm state
                sm{state_id_idx} = sm_par.state_ids(sm{state_idx}); % update sm state_id
                sm_trace = horzcat(sm_trace,[Trace(1,i)*T/N_trials;sm{state_idx}]); % final state
            end
        end
        % Stop in case of photobleaching
        if min(abs(sm{state_idx}-[sm_par.bleached_states]))==0
            sm{bleached_idx}=im_par.current_frame;
            sm{blinked_idx}=0;
            break;
        end
    end
    
    sm_trace = horzcat(sm_trace,[T ;sm{state_idx}]); % close up the trace

    if (sm{activated_idx}==0 && sm{state_idx}>=sm_par.converted_state); sm{activated_idx}=im_par.current_frame; end % update activation state if photoconversion occurred
    
else % no change in photophysical state
    sm_trace = horzcat(sm_trace,[T ;sm{state_idx}]); % simply close up trace
end

sm{state_trace_idx}=sm_trace;  % updated sm state trace
fluo_trace=sm_trace; % define fluorescent trace
fluo_trace(2,:)=sm_par.state_ids(fluo_trace(2,:)); % Replace by state_ids for rapidly exchanging states
fluo_trace(2,~ismember(fluo_trace(2,:),[sm_par.fluorescent_states]))=0;

sm{fluo_trace_idx}=fluo_trace;  % updated fluorescent state

%update blinking status at the end of the trace
if min(abs(sm_par.state_ids(sm_trace(2,end))-[sm_par.fluorescent_states]))~=0 && sm{bleached_idx}==0
    sm{blinked_idx}=1;
else
    sm{blinked_idx}=0;
end

%update total traces
sm{tot_state_trace_idx}=update_tot_trace_pct(sm([state_trace_idx,tot_state_trace_idx]), im_par);

%update Ns and sampling rate
if im_par.during_frametime==1
    sm{Ns_idx}(end,:)=Ns;
    sm{sampling_rate_idx}(1)=sampling_rate;
else
    sm{Ns_idx}(1,:)=Ns;
    sm{sampling_rate_idx}(2)=sampling_rate;
end

if plot_traces==1
    figure(1); clf
    if plot_fluo_trace==0
        p=plot(sm_trace(1,:),sm_trace(2,:));
    else
        p=plot(fluo_trace(1,:),fluo_trace(2,:));
    end
    
    set(p,'LineWidth',1.5);
    set(p,'Color','blue');
    ylim([0 max(ending_states(2,:))]);
    xlabel('Time [s]','fontsize',8,'fontweight','b')
    ylabel('State [AU]','fontsize',8,'fontweight','b')
    title(gca,'Single Molecule Trace','FontWeight','bold');
    disp(['Single molecule id:' , num2str(sm{id_idx})]);
    disp(['Photophysical state at end of trace:' , num2str(sm{state_idx})]);
    drawnow
    input('Ok ? ','s');
end






