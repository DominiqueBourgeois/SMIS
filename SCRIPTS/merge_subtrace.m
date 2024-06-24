function merged_trace=merge_subtrace(trace, subtrace)
%
% INPUTS:
%   trace: the current sm state trace
%   subtrace: the subtrace to append to trace
%	im_par: the imaging parameters

% MODIFICATION HISTORY:
%	D.Bourgeois, July 2020.

if ~isempty(trace)
    %     if trace(2,end)==subtrace(2,1) % If no change in state remove last point in trace
    merged_trace=horzcat(trace(:,1:end-1),subtrace);
    %     else %
    %         disp('h')
    %         merged_trace=horzcat(trace,subtrace);
    % end
else
    merged_trace=subtrace;
end




