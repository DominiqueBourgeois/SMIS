function sms=update_tot_fluo_trace(sms, sm_par)

% Retrieve the tot fluo traces from sm.tot_state_trace
% INPUTS:
%   sms: the single molecules
%	sm_par: the sm parameters
%	im_par: the imaging parameters
% OUTPUT
%   sms: the updated sms
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2020.

plot_trace=0; % Set to 1 to plot the trace

for i=1:numel(sms)
    for j=1:size(sms(i).sm,2)
        %Get the fluorescence trace
        tot_fluo_trace=sms(i).sm(j).tot_state_trace;
        tot_fluo_trace(2,:)=sm_par(i).state_ids(tot_fluo_trace(2,:)); % Replace by state_ids for rapidly exchanging states
        tot_fluo_trace(2,~ismember(tot_fluo_trace(2,:),[sm_par(i).fluorescent_states]))=0;
        sms(i).sm(j).tot_fluo_trace=tot_fluo_trace;
        
        if plot_trace==1
            t=tot_fluo_trace(1,:);
            s=tot_fluo_trace(2,:);
            %             tr2=[circshift(t,-1);s];
            %             tr2=tr2(:,1:end-1);
            %             tr3=horzcat(tr2,tot_fluo_trace);
            %             [~,s_order]=sort(tr3(1,:));
            %             tr3(:,:)=tr3(:,s_order(1:end));
            %             plot(tr3(1,:),tr3(2,:))
            plot(t,s)
            drawnow;
            xlabel('Time [s]');
            ylabel('Fluorescence state [s]');
            %             s3=tr3(2,:);
            %             ylim([0 1.1*max(s3)])
            ylim([0 1.1*max(s)])
        end
    end
end


% %OLD CODE
% bleach_frame=sms(i).sm(j).bleached; % Check if molecule was bleached and only consider trace until bleaching
% if bleach_frame > 0
%     all_state_traces=sms(i).sm(j).all_state_traces(:,1:bleach_frame);
%     N=bleach_frame;
% else
%     all_state_traces=sms(i).sm(j).all_state_traces;
%     N=im_par.current_frame;
% end
% % Calculate the start time for all traces
% if im_par.addtime>0 % For addtime
%     T_add=1e-3*(im_par.frametime+im_par.addtime)*(0:(N-1));
%     for k=1:N
%         all_state_traces{1,k}(1,:)=all_state_traces{1,k}(1,:)+T_add(k);
%     end
%     tot_trace_addtime=cell2mat(all_state_traces(1,:));
% end
% if im_par.frametime>0 % And for frametime
%     T_frame=1e-3*(im_par.frametime+im_par.addtime)*(0:(N-1))+1e-3*im_par.addtime;
%     for k=1:N-1 % Do it for 1:N-1, as the last may be empty
%         all_state_traces{end,k}(1,:)=all_state_traces{end,k}(1,:)+T_frame(k);
%     end
%     if ~isempty(all_state_traces{end,N})
%         all_state_traces{end,N}(1,:)=all_state_traces{end,N}(1,:)+T_frame(N);
%     end
%     tot_trace_frametime=cell2mat(all_state_traces(end,:));
% end
% if im_par.addtime>0 %Reorder in time if addtime
%     tot_trace=horzcat(tot_trace_addtime,tot_trace_frametime);
%     [~, si]=sort(tot_trace(1,:));
%     sms(i).sm(j).tot_state_trace=tot_trace(:,si);
% elseif im_par.frametime>0
%     sms(i).sm(j).tot_state_trace=horzcat(tot_trace_addtime,tot_trace_frametime);
% else
%     sms(i).sm(j).tot_state_trace=[;];
% end




