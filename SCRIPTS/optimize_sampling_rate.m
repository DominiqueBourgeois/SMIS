function [sampling_rate, trans_p] = optimize_sampling_rate(sampling_rate, sm_par, trans_p, minimum_oversampling, min_samples_per_frame, Ns, T)

%
% PURPOSE:
%   Optimize the sampling rate to minimize computing time
%
% INPUTS:
%   sampling_rate: [s-1] input sampling_rate
%	sm_par: the sm parameters
%   trans_p: transition probability matrix 
%   minimum_oversampling: minimum_oversampling
%   min_samples_per_frame: minimum # of samples per frame 
%   Ns: # of absorbed photons
%   T: time window to study in [s]

%
% OUTPUTS:
%   sampling_rate: [s-1] updated sampling_rate
%   trans_p: updated transition probability matrix 
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2022.

% It is also possibly to decrease the sampling rate if the
% option im_par.optimize_sampling_rate has been chosen. In that case, set
% the threshold a little bit lower (0.08). Run the while loop as above
optimized_sampling=0;
if (max(trans_p(:)) < 0.8/minimum_oversampling) && sampling_rate*T>=min_samples_per_frame
    while max(trans_p(:)) < 0.8/minimum_oversampling && sampling_rate*T>=min_samples_per_frame && optimized_sampling==0
        old_sampling_rate=sampling_rate;
        %new_sampling_rate=minimum_oversampling*old_sampling_rate*max(trans_p(:));
        new_sampling_rate=max([minimum_oversampling*old_sampling_rate*max(trans_p(:)),min_samples_per_frame/T]);
        % There must be an integer # of sampling in T
        required_number_of_samples=ceil(new_sampling_rate*T);
        new_sampling_rate=required_number_of_samples/T;
        if new_sampling_rate < old_sampling_rate
            sampling_rate=new_sampling_rate;
            %and redo the job
            trans_p_k=sm_par.trans_k/sampling_rate;
            Ns=Ns*old_sampling_rate/new_sampling_rate;
            trans_p_q = 1-bsxfun(@power,1-sm_par.trans_q,Ns'); % transition matrix for photo activated processes (each row should be powered by a corresponding N)
            trans_p=trans_p_k+trans_p_q-trans_p_k.*trans_p_q;
            trans_p=trans_p'; % Transpose matrix (for further processing with get_state_evolution.m)
        else
            optimized_sampling=1; % Cannot go further
        end
    end
end
