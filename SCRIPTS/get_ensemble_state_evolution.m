function [sp, sm_par]=get_ensemble_state_evolution(sp, lasers, sm_par, im_par)

%
% PURPOSE:
%   Get the evolution of photophysical state at the ensemble level
%
% INPUTS:
%   sp: the sub population (with coordinates on high-resolution image) in raster units
%   lasers: the lasers
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%   w_events: structure array of indices were transitions might be detected (only the variable is iput, its content is derived in the present routine)
%   starting_states: array where starting states for the different transitions are specified
%   ending_states: array where ending states for the different transitions are specified
%
% OUTPUTS:
%   sp: the sub population updated for photophysical state
%	sm_par: the sm parameters eventually updated for sampling rate
%   sm_par.sampling_rate(1) is sampling_rate for frametime and sm_par.sampling_rate(2) is sampling_rate for addtime
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2021.


% plot_traces=0; % Set to 1 if every trace to be displayed
% plot_fluo_trace=0; % Set to 1 to see fluorescent state appears

% get general parameters
% Time window to study in [s]
if im_par.during_frametime==1
    T=1e-3*im_par.frametime;
    %Depending wether used_sampling_rate has already been calculated choose
    %initial sampling rate
    if ~isnan(sp.sampling_rate(1))
        sampling_rate=sp.sampling_rate(1);
    else
        sampling_rate=sm_par.start_sampling_rate(1);
    end
else
    T=1e-3*im_par.addtime;
    if ~isnan(sp.sampling_rate(2))
        sampling_rate=sp.sampling_rate(2);
    else
        sampling_rate=sm_par.start_sampling_rate(2);
    end
end

minimum_oversampling=im_par.minimum_oversampling;

% Trace=[;]; % initialize trace

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

% So we now calculate the transition rate matrix

% transition rate matrix for thermal processes, rates given in [s-1]
% trans_p_k=sm_par.trans_k/sampling_rate;
trans_k=sm_par.trans_k;

% transition rate matrix for photo activated processes
% # of absorbed photons per sampling time in each starting photoactive state
Ns=get_number_of_absorbed_photons_ensemble(sp, lasers, sm_par, im_par);
% equivalent probability transition matrix, rates given in [s-1], see discussion above
% trans_p_q = 1-bsxfun(@power,1-sm_par.trans_q,Ns'); % transition matrix for photo activated processes (each row should be powered by a corresponding N)
trans_q = bsxfun(@times,sm_par.trans_q,Ns'); % rate matrix for photo activated processes 


% overall probability tansition matrix, rates given in [s-1]
% trans_p=trans_p_k+trans_p_q-trans_p_k.*trans_p_q;
rate_matrix=trans_k+trans_q;

% to gain speed we should now define which state can be accessed from what
% state
% [trans_p, allowed_n_trans, allowed_trans_indices]=get_allowed_transitions(trans_p, sp, sm_par);

% Check sampling rate. Sampling rate should be at least minimum_oversampling times more than actual rates.
% Actual rates [s-1] = trans_p*sampling_rate so that trans_p*sampling_rate must be < sampling rate/minimum_oversampling
% This translates into max(trans_p) < 1/minimum_oversampling. If not increase sampling_rate

if max(rate_matrix(:)) > sampling_rate/minimum_oversampling
    % Run the following as many times as necessary
    while max(rate_matrix(:)) > sampling_rate/minimum_oversampling
        new_sampling_rate=minimum_oversampling*max(rate_matrix(:));
        
        % There must be an integer # of sampling in time window T
        required_number_of_samples=ceil(new_sampling_rate*T);
        new_sampling_rate=required_number_of_samples/T;
        sampling_rate=new_sampling_rate;
    end
%     optimized_sampling=1; % Set this flag for next optimization option below
    if im_par.during_frametime==1
        disp(['Sampling rate during frametime has to be increased to: ',num2str(sampling_rate),' Hz'])
    else
        disp(['Sampling rate during addtime has to be increased to: ',num2str(sampling_rate),' Hz'])
    end
else
%     optimized_sampling=0; % Set this flag for next optimization option below
end

% Now do the calculation
n=sampling_rate*T; % Number of steps
dt=1/sampling_rate; % Time step
S=[sp.p]'; % set the current populations
for k=1:n
    S= S+dt*(-sum(rate_matrix,2).*S+rate_matrix'*S);
    
    % Equilibrate states in fast equilibrium
    if sm_par.pH_sensitivity==1
       for i=1:sm_par.n_fluorescent_states
           tmp_C=S(sm_par.fluorescent_states(i))+S(sm_par.associated_dark_states(i));
           S(sm_par.fluorescent_states(i))=sm_par.fluorescent_fraction(i)*tmp_C;
           S(sm_par.associated_dark_states(i))=(1-sm_par.fluorescent_fraction(i))*tmp_C;
       end
    end    
    
    if im_par.during_frametime==1
        sp.det_p(:,im_par.current_frame)=sp.det_p(:,im_par.current_frame)+S; % Integrate the population over time
    end
end
sp.p=S'; % Assign to current populations

if im_par.during_frametime==1
    sp.det_p(:,im_par.current_frame)=sp.det_p(:,im_par.current_frame)/n; % Divide by the # of steps

    %Eventually reassign states in rapid equilibrium
    if sm_par.pH_sensitivity==1
       for i=1:sm_par.n_fluorescent_states
           sp.det_p(sm_par.fluorescent_states(i),im_par.current_frame)=...
               sp.det_p(sm_par.fluorescent_states(i),im_par.current_frame)+sp.det_p(sm_par.associated_dark_states(i),im_par.current_frame);
           sp.det_p(sm_par.associated_dark_states(i),im_par.current_frame)=0;
       end
    end    


end

%update Ns and sampling rate
if im_par.during_frametime==1
    sp.Ns(end,:)=Ns;
    sp.sampling_rate(1)=sampling_rate;
else
    sp.Ns(1,:)=Ns;
    sp.sampling_rate(2)=sampling_rate;
end







