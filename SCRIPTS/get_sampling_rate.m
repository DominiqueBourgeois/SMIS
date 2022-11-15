function sm_par=get_sampling_rate(sm, lasers, sm_par, im_par)

%
% PURPOSE:
%   Get initial sampling rate at the average molecule position 
%
% INPUTS:
%   sm: the average single molecule (with coordinates on high-resolution image) in raster units
%   lasers: the lasers
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%	sm_par: the sm parameters eventually updated for sampling rate
%   sm_par.start_sampling_rate(1) is sampling_rate for frametime and sm_par.start_sampling_rate(2) is sampling_rate for addtime
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2022, adapted to parallel computing

% Time window to study in [s]
if im_par.during_frametime==1
    T=1e-3*im_par.frametime;
    sampling_rate=sm_par.sampling_rate(1);
else
    T=1e-3*im_par.addtime;
    sampling_rate=sm_par.sampling_rate(2);
end

minimum_oversampling=im_par.minimum_oversampling;
min_samples_per_frame=im_par.min_sampling_points;

% *** Compute probability of transitions *** (see get_state_evolution)
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
Ns=get_number_of_absorbed_photons(sm, lasers, sm_par, im_par, sampling_rate);

% equivalent probability transition matrix, rates given in [s-1], see discussion above
trans_p_q = 1-bsxfun(@power,1-sm_par.trans_q,Ns'); % transition matrix for photo activated processes (each row should be powered by a corresponding N)

% overall probability tansition matrix, rates given in [s-1]
trans_p=trans_p_k+trans_p_q-trans_p_k.*trans_p_q;

% to gain speed we should now define which state can be accessed from what
% state
[trans_p, ~, ~]=get_allowed_transitions(trans_p, sm, sm_par);

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
        [trans_p, ~, ~]=get_allowed_transitions(trans_p, sm, sm_par);
    end
    optimized_sampling=1; % Set this flag for next optimization option below
    if im_par.during_frametime==1
        disp(['Sampling rate during frametime has to be increased to: ',num2str(sampling_rate),' Hz'])
    else
        disp(['Sampling rate during addtime has to be increased to: ',num2str(sampling_rate),' Hz'])
    end
else
    optimized_sampling=0; % Set this flag for next optimization option below
end

% Eventually optimize sampling rate
if im_par.optimize_sampling_rate==1 && optimized_sampling==0
    [sampling_rate, ~] = optimize_sampling_rate(sampling_rate, sm_par, trans_p, minimum_oversampling, min_samples_per_frame, Ns, T);
end

%update sampling rate
if im_par.during_frametime==1
    sm_par.start_sampling_rate(1)=sampling_rate;
    disp(['Sampling rate for fluorophore ',sm_par.fluorophore_name,' during frametime set to: ',num2str(sampling_rate),' Hz'])
else
    sm_par.start_sampling_rate(2)=sampling_rate;
    disp(['Sampling rate for fluorophore ',sm_par.fluorophore_name,' during addtime set to: ',num2str(sampling_rate),' Hz'])
end

