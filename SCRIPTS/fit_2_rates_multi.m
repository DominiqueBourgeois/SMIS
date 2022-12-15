function y = fit_2_rates_multi(a, x)
%
% MODIFICATION HISTORY:
%       D.Bourgeois, February 2012.
%-
global start start_common w_variables w_variables_common nsets datasets_lengths;

a2 = start;
a2_common = start_common;

if ~isempty(w_variables_common)
    a2_common(w_variables_common) = a(1:length(w_variables_common));
end
if ~isempty(w_variables)
    a2(w_variables) = a(length(w_variables_common)+1:end);
end

s_a2=size(a2,1);
k=1;

for i=1:nsets
    loc_x=x(k:k+datasets_lengths(i)-1);
    % We have to define the 5 param [a0_1(i), a0_2(i), offset(i), rate1, rate2];
    loc_a2=[a2(1:s_a2,i)',a2_common'];
    loc_y = get_2_rates_model(loc_a2,loc_x);
    if i==1 ; y=loc_y; end
    if i>1
        y=vertcat(y,loc_y);  %#ok<AGROW>
    end
    k=k+datasets_lengths(i);
end
end



