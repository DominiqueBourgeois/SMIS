function y = fit_exp_function_multi_aniso_linear(a, x)
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

k=1;

for i=1:nsets
    loc_x=x(k:k+datasets_lengths(i)-1);
    % We have to define the 4 param [a0(i), k_a(i), offset(i), phi, beta];
    loc_a2=[a2(1:size(a2,1),i)',a2_common'];
    loc_y = get_exp_aniso_linear_model(loc_a2,loc_x);
    if i==1 ; y=loc_y; end;
    if i>1 ; y=vertcat(y,loc_y); end;
    k=k+datasets_lengths(i);
end
end



