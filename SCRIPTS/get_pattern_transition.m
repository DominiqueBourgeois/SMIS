function [S2, p_id] = get_pattern_transition(P,V, D_ras,S,K,C,dt, sm_par, im_par)

% PURPOSE:
% Get new diffusion state for a molecule initially in diffusion state S
% with state transition matrix K containing rates [s-1]
%
% INPUTS:
%   P: Current XYZ coordinate of SM before current diffusion + velocity step
%   V: Current Velocity
%   D_ras: Current diffusion coefficients in raster^.s-1
%   S: Starting diffusion state
%   K(S1,S2) is the exchange rate between state S1 and S2
%   C: confinement sub-pattern ids for diffusion states
%   dt: time interval [s]
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%	S2: new diffusion state
%	p_id: new sub_pattern id
%
% METHOD:
% Relation between rate [s-1) and probability of transition
% k=p*n with k = rate [s-], p = probability (or yield) and n = # of trials per second
% so p = k/n
% Exemple: rate = 10s-1 and we test every dt [sec]. If we test for T [sec],
% we make T/dt attempts which a chance of success each time of k/(1/dt) = k*dt
% so we test N=T/dt random numbers between 0 and 1 and look those that are <= k*dt
% So here if we test for N=100 (T=frametime, dt=frametime/N), we check if 100 random number are <=
% k*frametime/N
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15
% Modified on November 8 2019 (bug)
%	D.Bourgeois, November 2019: version > simulate_palm_vsn15.3
%	D.Bourgeois, December 2019: version > simulate_palm_vsn15.3
%	D.Bourgeois, December 2020: corrected bug for indexing patterns:
%	sm_par.w_patterns(C_ES+1).w replaced by sm_par.w_patterns(ES(i)).w as
%	C_ES is the subpatern id (ie: 255), not the id of the subpattern (ie
%	between 1 and 3 if 3 subpatterns)
%	D.Bourgeois, November 2020: version > simulate_palm_vsn16.3 (introduce
%	Velocity)

if isempty(K); S2=S; p_id=C(S); return; end % No change in state in that case

N=100; % We will make N tests per molecule per frametime

Starting_States=K(:,1); % these are the starting states
Ending_States=K(:,2); % These are the ending states
Rates=K(:,3)'; % These are the rates, transpose the Matrix
w=find(Starting_States==S); % Look which transitions deal with our current starting state
R=Rates(w);

Trans_indices=find(R>0); % Indices of possible transitions
ES=Ending_States(w(Trans_indices)); % These are the possible actual ending states
N_trans=numel(Trans_indices); % Number of possible transitions
w_event=N+1; % define something > N
Test=rand(N_trans,N); % series of random numbers
S2=[];

%Handle Diff independant transitions
if sm_par.DIT(S)==1
    force_change=1;
else
    force_change=0;
end


for i=1:N_trans
    w = find(Test(i,:) <= dt/N*R(Trans_indices(i)),1); % We just look for the first event here to make it simpler
    if ~isempty(w) % Possible change in diffusion identified
        %         if w<w_event || w==1 % The case w=1 should be considered if several highly probable transitions occur at w=1.
        if w<=w_event % The case w=w_event should be considered if several transitions occur at w.
            % Then the choice between involved transitions will be made at the end
            % Examine if this change is possible
            C_ES=C(ES(i)); % Subpattern id of ending state
            if C(S)==C_ES  % No change in confinement, so no limitation here
                possible_change_in_S=1;
            elseif force_change==1 % Case of diffusion independant transition
                possible_change_in_S=1;
            else
                if im_par.simul_3D==0 % 2D case
                    % Get radius of exploration for current diffusion state
                    % D_ras in pix^2/s, ie D_ras*dt is the explored
                    % surface in dt, ie a disk of radius pi*r^2 =
                    % D_ras*dt
                    r=sqrt(D_ras*dt/pi);
                elseif im_par.simul_3D==1 % 3D case
                    % D_ras*frametime is the explored
                    % volume in frametime, ie a sphere of radius 4/3*pi*r^3 =
                    % D_ras*frametime
                    r=(3/(4*pi)*D_ras*dt).^(1/3);
                end
                np_id=sm_par.n_sp_id==sm_par.D_confined(ES(i)); % Pattern id for new diffusion state
                possible_change_in_S = search_close_subpatterns(P,V*dt,sm_par.w_patterns(np_id).w,r,im_par); % Changed 17/01/21
            end
            
            % If the change is possible, go for it !
            if possible_change_in_S==1
                %                 w_event=w; % If more transitions are to be tested, they should occur before the current successful transition, so we decrease w_event to the value of w
                %                 if w==1 % This might happen if several transitions are allowed with w=1
                if w==w_event % If several transitions occur at w
                    S2=[S2,ES(i)];  %#ok<AGROW>
                else
                    S2=ES(i);
                end
                w_event=w; % If more transitions are to be tested, they should occur before the current successful transition, so we decrease w_event to the value of w
            end
        end
    end
end

if isempty(S2)
    S2=S;
elseif numel(S2)>1
    S2=S2(randi([1,numel(S2)])); % If several transitions occur concomitantly, pick randomly one of them
end
p_id=C(S2); % Update pattern_id for SM

end

