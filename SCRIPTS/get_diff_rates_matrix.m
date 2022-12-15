function K=get_diff_rates_matrix(N,diff_exchange_rates)
% Get state transition matrix K containing rates [s-1]
% K(S1,S2) is the exchange rate between state S1 and S2
% N is the # of diffusion states
% diff_exchange_rates is the exchange rates between diffusion coefficients: array of [3 X Nx(N-1)] with N=(length(diff_coefficients)]
% Each series of 3 numbers consists of "starting state", "ending state", "exchange rate [s-1] between starting and ending states

s=size(diff_exchange_rates,1);
if s>0
    K=zeros(N,N) ; %  transition Probabilité matrice
    for k=1:s
        K(diff_exchange_rates(k,1),diff_exchange_rates(k,2))=diff_exchange_rates(k,3);
    end
else
    K=[];
end
end
