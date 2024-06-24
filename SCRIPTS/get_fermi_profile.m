function [fermi_profile, frame_series] = get_fermi_profile(n_frames, n_frames_fermi, frame_falloff, plot_figure)

%+
%
% MODIFICATION HISTORY:
%       D.Bourgeois, Décembre 2014 .
%-
% Get Fermi profile for 405 nm laser activation according to Lee et al,
% PNAS 2012

frame_series=1:n_frames;

T=frame_falloff;
tF=n_frames_fermi;
t=frame_series;
f=exp(-(t-tF)/T);
fermi_profile= f./((1+f).*(log(1+f)));

if max(fermi_profile)>1
    plat_start = find(fermi_profile > 0.9999, 1, 'first');
    fermi_profile(plat_start:end)=1;
end

if plot_figure==1
    pf=100./(1+f);
    figure;
    hold on
    title('Fermi plot: ','FontWeight','bold');  
    subplot(2,1,1)
    plot(frame_series,fermi_profile);
    xlabel('Frame series','fontsize',8,'fontweight','b')
    ylabel('Fermi profile','fontsize',8,'fontweight','b')
    subplot(2,1,2)
    plot(frame_series,pf,'color','r');
    xlabel('Frame series','fontsize',8,'fontweight','b')
    ylabel('PDF','fontsize',8,'fontweight','b')
%     legend('Fermi profile','PDF [Red]');
    disp(['Fraction of power used in first frame [%]: ', num2str(100*fermi_profile(1))]);
    disp(['Fraction of power used in last frame [%]: ', num2str(100*fermi_profile(end))]);
    hold off
end


