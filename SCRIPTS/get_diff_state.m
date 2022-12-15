function sm=get_diff_state(sm, K, confined)

% PURPOSE:
% Get initial diffusion state for a dye
%
% INPUTS:
%	sm: the single molecules
%	K: the diffusion coefficient rate matrix 
%   confined: the pattern id of the diffusion states 
%
% OUTPUTS:
%	sm: the updated single molecules
%
% MODIFICATION HISTORY:
%	D.Bourgeois, December 2019: version > simulate_palm_vsn15.3
%   D.Bourgeois, January 2022: Make sure there are no leftover molecules with no assigned diffusion state
%   D.Bourgeois, February 2022: Corrected bug: should be
%   d_ini=equilibrate_initial_values(K) (not K') !!


% Go through all confinement states
for k=0:max(confined)
    D_state_ids=find(confined==k);% This gives us the indices of diffusion states that are confined in subpattern k
    if ~isempty(D_state_ids) % Then do something
        if length(D_state_ids)==1 % Just 1 diffusion state in subpattern
            [sm([sm.c_sp]==k).diff_state]=deal(D_state_ids);
        else
            % In this case there are several diffusion states in
            % the sub pattern
            % The best that we can do is to distribute the
            % molecules according to the kinetic scheme, But
            % ignoring spatial distribution of the patterns.
            d_ini=equilibrate_initial_values(K); % Equilibrate inside that subpattern.
            if isempty(d_ini) % Case of no exchange
                d_ini=ones(1,numel(confined));
            end
            d_ini=d_ini(D_state_ids)/sum(d_ini(D_state_ids)); % This gives the fraction of molecules with the associated diffusion states
            w=find([sm.c_sp]==k);
            if ~isempty(w)
                n_mol=numel(w);
                r=randperm(n_mol,n_mol); % reorder molecules randomly
                j0=1; % indices
                for j=1:length(D_state_ids)
                    n=round(n_mol*d_ini(j));
                    [sm(w(r(j0:min([j0+n-1, n_mol])))).diff_state]=deal(D_state_ids(j));
                    j0=j0+n;
                end
            end
            %Make sure there are no leftover molecules with no assigned diffusion state
            emptyIndex = find(arrayfun(@(MyStruct) isempty(MyStruct.diff_state),sm(w)));
            if ~isempty(emptyIndex)
                max_pos=find(d_ini==max(d_ini),1,'first');
                [sm(w(emptyIndex)).diff_state]=deal(D_state_ids(max_pos));
            end
        end
    end
end
end
