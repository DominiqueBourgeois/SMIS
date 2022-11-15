function y = fit_N_blink_hist(a, x)
%
% MODIFICATION HISTORY:
%       D.Bourgeois, July 2013.
%-
global start w_variables

a2 = start;

if ~isempty(w_variables)
    a2(w_variables) = a(1:length(w_variables));
end

y = get_N_blink_model(a2,x);
end



