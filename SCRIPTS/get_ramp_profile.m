function [ramp_profile, frame_series] = get_ramp_profile(n_frames, min_power, max_power, step_size, plot_figure)

%+
%
% MODIFICATION HISTORY:
%       D.Bourgeois, Décembre 2014 .
%-
% Get Ramp profile for 405 nm laser activation

frame_series=1:n_frames;
ramp_profile=zeros(1,n_frames);
% n_steps=fix(n_frames/step_size);
p=min_power;
s=0;
step_size=round(step_size);
for k=(n_frames-step_size):n_frames
    if mod(k,step_size)==1
        s=k;
    end
end

if step_size==1
    ramp_profile=min_power+(0:(n_frames-1))/(n_frames-1)*(max_power-min_power);
else
    
    for k=1:n_frames
        if mod(k,step_size)==1
            %    p=min_power+(k-1)*step_size/n_steps*(max_power-min_power);
            %    p=min_power+(k-1)/(n_frames-step_size)*(max_power-min_power);
            p=min_power+(k-1)/(s-1)*(max_power-min_power);
        end
        ramp_profile(k)=p;
    end
end


if plot_figure==1
    figure;
    title('Ramp plot: ','FontWeight','bold');
    plot(frame_series,ramp_profile);
    xlabel('Frame series','fontsize',8,'fontweight','b')
    ylabel('Ramp profile','fontsize',8,'fontweight','b')
    disp(['Power used in first frame [mW]: ', num2str(ramp_profile(1))]);
    disp(['Power used in last frame [mW]: ', num2str(ramp_profile(end))]);
    hold off
end


