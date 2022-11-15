function display_par=plot_final_photophysics(im_par, sm_par, display_par)

% PURPOSE:
%	Plot final statistics on photophysical state of molecules
%
% INPUTS:
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   h: the figure id
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019: version > simulate_palm_vsn15

frame=im_par.current_frame;
if frame>1
    
    % Create a figure if necessary
    create_figure=1;
    % The following code ensures that the right figure is addressed
    % Indeed, it might be that the figure has been created before, then deleted, and the same figure number has been re-created for another purpose and is still open
    
    % Take a look at the current figures open
    figHandles = findobj('Type', 'figure'); % Get the array of all figures open
    n_fig=size(figHandles,1);
    for k=1:n_fig
        if strcmp(figHandles(k).Name,'Photophysics Summary')==1 % Our figure already exists
            create_figure=0;
            display_par.photophysics_figure_number=figHandles(k).Number;
        end
    end
    
    
    % create a figure if necessary
    if create_figure==1
        w_xsize=display_par.scrsz(3)*3/5;
        w_ysize=display_par.scrsz(4)*3/5;
        photophysics_figure=figure;
        display_par.photophysics_figure_number=photophysics_figure.Number;
        colormap(gray);
        set(photophysics_figure,'Position',[10 display_par.scrsz(4)/10 w_xsize w_ysize])
        set(photophysics_figure,'color','white');
        set(photophysics_figure,'Name','Photophysics Summary'); 
        set(photophysics_figure,'NumberTitle','off');
    else
        figure(display_par.photophysics_figure_number);
        clf;
    end
    
    
    n_fluorophores=size(sm_par,2);
    
    
    x=1:frame;
    
    for k=1:n_fluorophores
        photophysical_status=sm_par(k).photophysical_status;
        % Bleaching
        subplot(n_fluorophores,3,1+3*(k-1));
        
        hold on;
        y1=photophysical_status.n_bleached(1:frame);
        y2=photophysical_status.cum_n_bleached(1:frame);
        p=plot(x,y1);
        set(p,'LineWidth',1.5);
        set(p,'Color','blue');
        line(x,y2,'color','g','LineWidth',1.5);
        ylim([0 max([1.2*max([max(y1) max(y2)]) 1])]);
        xlim([1 frame]);
        xlabel('Frame #','fontsize',8,'fontweight','b')
        ylabel('# of molecules','fontsize',8,'fontweight','b')
        title(gca,['Bleaching: ', sm_par(k).fluorophore_name],'FontWeight','bold');
        legend('Current','Cumulated');
        
        % Activation
        subplot(n_fluorophores,3,2+3*(k-1));
        hold on;
        y1=photophysical_status.n_activated(1:frame);
        y2=photophysical_status.cum_n_activated(1:frame);
        p=plot(x,y1);
        set(p,'LineWidth',1.5);
        set(p,'Color','blue');
        line(x,y2,'color','g','LineWidth',1.5);
        ylim([0 max([1.2*max([max(y1) max(y2)]) 1])]);
        xlim([1 frame]);
        xlabel('Frame #','fontsize',8,'fontweight','b')
        ylabel('# of molecules','fontsize',8,'fontweight','b')
        title(gca,['Activation: ', sm_par(k).fluorophore_name],'FontWeight','bold');
        legend('Current','Cumulated');
        
        % Blinking
        subplot(n_fluorophores,3,3+3*(k-1));
        hold on;
        y1=photophysical_status.n_blinked(1:frame);
        p=plot(x,y1);
        set(p,'LineWidth',1.5);
        set(p,'Color','blue');
        ylim([0 max([1.2*max(y1) 1])]);
        xlim([1 frame]);
        xlabel('Frame #','fontsize',8,'fontweight','b')
        ylabel('# of molecules','fontsize',8,'fontweight','b')
        title(gca,['Blinking: ', sm_par(k).fluorophore_name],'FontWeight','bold');
        legend('Current');
    end
end