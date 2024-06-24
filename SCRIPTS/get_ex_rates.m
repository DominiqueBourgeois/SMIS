function K=get_ex_rates(exchange_matrix)

% Inverse routine as compared to get_diff_rates_matrix.m
% Get array of [3 X Nx(N-1)] with N=(length(diff_coefficients)]
% from state transition matrix K containing rates [s-1]
% exchange_matrix(S1,S2) is the exchange rate between state S1 and S2
% K is the exchange rates between diffusion coefficients: array of [3 X Nx(N-1)] with N=(length(diff_coefficients)]
% Each series of 3 numbers consists of "starting state", "ending state", "exchange rate [s-1] between starting and ending states

if ~isempty(exchange_matrix)   
    [w_row,w_col]=find(exchange_matrix>0);   
    n_transitions=numel(w_col);
    if n_transitions>0
            K=zeros(n_transitions,3);
            for i=1:n_transitions
                K(i,:)=[w_row(i),w_col(i),exchange_matrix(w_row(i),w_col(i))] ;
            end
    else
        K=[];
    end
else
    K=[];
end
