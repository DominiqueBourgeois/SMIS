function [trans_p, allowed_n_trans, allowed_trans_indices]=get_allowed_transitions(trans_p, sm, sm_par)

% PURPOSE:
%   Depending on the state of a sm, return the number and the indices of
%   the allowed state transitions.
%   if the molecule is not photoconverted: everything can happen
%   if the molecule is photoconverted: only those states corresponding to a
%   photoconverted state can occur
%
% INPUTS:
%   sm: the single molecule 
%	sm_par: the sm parameters
%	trans_p: the transition rate matrix: trans_p(i,j) is the transition
%	from state i to state j
%
% OUTPUTS:
%	trans_p: the updated and transposed transition rate matrix
%   allowed_n_trans: the number of allowed transition
%   allowed_trans_indices: the indices of the allowed transitions
%
% MODIFICATION HISTORY:
%	D.Bourgeois, July 2019.

%test if sm is activated
if sm.activated>0
    trans_p(1:sm_par.converted_state-1,:)=0;
else
end
trans_p=trans_p';
allowed_trans_indices=find(trans_p>0);
allowed_n_trans=numel(allowed_trans_indices);
