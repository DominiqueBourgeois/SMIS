function [subTrace,sampling_rate]=get_fret_subtrace_pct(state,sm_par,N,sampling_rate,start_time,end_time, im_par, w_events, state_link)

%
% PURPOSE:
%   Get the first evolution of photophysical state of a molecule in a
%   defined initial state
%
% INPUTS:
%   state: current state of single molecule 
%	sm_par: the sm parameters
%   N: the # of photons per sampling time of the sm in state 
%   sampling_rate: the sm sampling rate
%   start_time: starting time of the subtrace
%   end_time: maximum ending time of the subtrace
%   im_par: imaging parameters
%   w_events: structure array of indices were transitions might be detected (only the variable is iput, its content is derived in the present routine)
%   state_link: structure array equivalent to starting_states and ending_states, but for submatrix of trans_k and trans_q with defined starting states  
%
% OUTPUTS:
%   subTrace: the computed subTrace
%   sampling_rate: the updated sampling rate
%
% MODIFICATION HISTORY:
%	D.Bourgeois, July 2019.
%	D.Bourgeois, July 2020.
%	D.Bourgeois, September 2022, adapted to parallel computing


plot_traces=0; % Set to 1 if every trace to be displayed
plot_fluo_trace=0; % Set to 1 to see fluorescent state appears
print_message=0; % Set to 1 to see message upon sampling state change

minimum_oversampling=im_par.minimum_oversampling;
min_samples_per_frame=im_par.min_sampling_points;

T_all=end_time-start_time; % Max duration of trace

Trace=[;]; % initialize subtrace

if im_par.during_frametime==1
    T=1e-3*im_par.frametime;
else
    T=1e-3*im_par.addtime;
end


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
trans_p_k=sm_par.trans_k(sm_par.state_ids==state,:)/sampling_rate;

% probability transition matrix for photo activated processes
% # of absorbed photons per sampling time in each starting photoactive state

% equivalent probability transition matrix, rates given in [s-1], see discussion above
trans_p_q = 1-bsxfun(@power,1-sm_par.trans_q(sm_par.state_ids==state,:),N'); % transition matrix for photo activated processes (each row should be powered by a corresponding N)

% overall probability tansition matrix, rates given in [s-1]
trans_p=trans_p_k+trans_p_q-trans_p_k.*trans_p_q;

% to gain speed we should now define which state can be accessed from what
% state
trans_p=trans_p';
allowed_trans_indices=find(trans_p>0);
allowed_n_trans=numel(allowed_trans_indices);

% Check sampling rate. Sampling rate should be at least minimum_oversampling times more than actual rates.
% Actual rates [s-1] = trans_p*sampling_rate so that trans_p*sampling_rate must be < sampling rate/minimum_oversampling
% This translates into max(trans_p) < 1/minimum_oversampling. If not increase sampling_rate

if max(trans_p(:)) > 1/minimum_oversampling
    % Run the following as many times as necessary
    while max(trans_p(:)) > 1/minimum_oversampling
        old_sampling_rate=sampling_rate;
        new_sampling_rate=minimum_oversampling*old_sampling_rate*max(trans_p(:));
        
        % There must be an integer # of sampling in time window T
        % Note: we use here the full frametime or addtime for T, even if we
        % process subtraces. We do not want the sampling rate to increase
        % when the length of the subtrace to study will shrink
        required_number_of_samples=ceil(new_sampling_rate*T);
        new_sampling_rate=required_number_of_samples/T;
        sampling_rate=new_sampling_rate;
        %and redo the job
        trans_p_k=sm_par.trans_k(sm_par.state_ids==state,:)/sampling_rate;
        N=N*old_sampling_rate/new_sampling_rate;
        
        trans_p_q = 1-bsxfun(@power,1-sm_par.trans_q(sm_par.state_ids==state,:),N'); % transition matrix for photo activated processes (each row should be powered by a corresponding N)
        trans_p=trans_p_k+trans_p_q-trans_p_k.*trans_p_q;
        
        trans_p=trans_p';
        allowed_trans_indices=find(trans_p>0);
        allowed_n_trans=numel(allowed_trans_indices);
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
    [sampling_rate, trans_p] = optimize_sampling_rate(sampling_rate, sm_par, trans_p, minimum_oversampling, min_samples_per_frame, N, T);
end


% prevent transition probabilities of 1, if they happen (singular case, see
% above)
trans_p(trans_p==1)=1-1e-5;

N_trials=round(T_all*sampling_rate); % Number of trials: we will do a test every 1/sampling_rate [s]

Test=rand(allowed_n_trans,N_trials); % series of random numbers
for i=1:allowed_n_trans
    w_events(i).w = find(Test(i,:) <= trans_p(allowed_trans_indices(i)));
    Trace=horzcat(Trace,[w_events(i).w;allowed_trans_indices(i)*ones(1,length([w_events(i).w]))]);  % concatenate all these events
    % for each event of the trace, associate index equal to Trans_indices(i)
end

subTrace=[start_time;state]; % initialize subtrace

% reorganize  trace in chronological order
if ~isempty(Trace)
    [~,ix]=sort(Trace,2);
    ix1=ix(1,:);
    Trace(:,:)=Trace(:,ix1(:));    % trace to process
        
    %process the 1st event that involves current state or states in
    %equilibrium with state
    S=state_link(state).starting_states(state_link(state).index==Trace(2,1)); % What is the starting state for this transition
    F=state_link(state).ending_states(state_link(state).index==Trace(2,1)); % What is the ending state for this transition
          
    subTrace = horzcat(subTrace,[start_time+Trace(1,1)*T_all/N_trials;S]); % initial state 
    subTrace = horzcat(subTrace,[start_time+Trace(1,1)*T_all/N_trials;F]); % final state  
    
else % no change in photophysical state
    subTrace = horzcat(subTrace,[end_time;state]); % simply close up trace
end



if plot_traces==1
    figure(1)
    if plot_fluo_trace==0
        p=plot(subTrace(1,:),subTrace(2,:));
    end
    set(p,'LineWidth',1.5);
    set(p,'Color','blue');
    ylim([0 max(state_link(state).ending_states)]);
    xlabel('Time [s]','fontsize',8,'fontweight','b')
    ylabel('State [AU]','fontsize',8,'fontweight','b')
    title(gca,'Single Molecule Trace','FontWeight','bold');
    disp(['Subtrace #:' , num2str(subtrace_id)]);
    disp(['New photophysical state:' , num2str(state)]);
    input('Ok ? ','s');
end






