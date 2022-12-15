function [mol_state, mol_state_trace]=check_mol_state(sm,sm_par, im_par, frame)
% NAME:
%	check_blinking_state
% PURPOSE:
%   Check the blinking state from a sm.
%
% INPUTS:
%
%
% OUTPUTS:
%
% MODIFICATION HISTORY:
%	D.Bourgeois, January 2015.

debug=0; % set to 1 to see messages
f=sm_par.fraction_state1;

% Treat the cases where all the population is in state 1 or state 2
if f==1
    mol_state=1;
    mol_state_trace=[];
    return
end

if f==0
    mol_state=2;
    mol_state_trace=[];
    return
end

if sm_par.exchange_rate==0
    disp('Error: The exchange rate cannot be 0 if the fraction of molecules in state 1 or 2 is not 0 or 1 !');
    return;
end

mol_state=sm.mol_state; % Initialize molecular state
state_exchange=0;
n_time_steps=(im_par.frametime+im_par.addtime)*1e+2; % number of time steps during entire frame in 10 us units
if sm.mol_state==1
    exchange_yield=sm_par.exchange_rate*1e-5; % Chance of exchange within time step
    test_exchange = rand(1, round(n_time_steps));
    w_exchange = find(test_exchange <= exchange_yield, 1);
    if ~isempty(w_exchange)
        mol_state=2;
        state_exchange=1;
    end
end
if sm.mol_state==2
    exchange_rate21=sm_par.exchange_rate*f/(1-f);
    if exchange_rate21>0
        exchange_yield=exchange_rate21*1e-5; % Chance of exchange within time step
        test_exchange = rand(1, round(n_time_steps));
        w_exchange = find(test_exchange <= exchange_yield, 1);
        if ~isempty(w_exchange)
            mol_state=1;
            state_exchange=1;
        end
    end
end

if state_exchange==1; % Update history of molecular state
    mol_state_trace=horzcat(sm.mol_state_trace,[frame;mol_state]);
else
    mol_state_trace=sm.mol_state_trace;
end
if debug==1
    disp(['Molecule was in molecular state: ', num2str(sm.mol_state)]);
    disp(['Molecule is now in molecular state: ', num2str(mol_state)]);
end




