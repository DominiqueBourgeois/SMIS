function plot_ensemble_data(ens, im_par, sm_par, smis_title)

% PURPOSE:
%	Plot final ensemble data
%
% INPUTS:
%	ens: the ensemble data
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%   title: title of SMIS simulation
%
% OUTPUTS:
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2021

n_fluorophores=numel(sm_par);
n_images=im_par.n_images;
append_figure=0; % Set to 1 to create new figures

% Pick the right figure
if append_figure==1
    if isempty(findall(0,'type','figure'))
        fig_num=1;
    else
        curr_fig=gcf;
        fig_num=curr_fig.Number+1;
    end
else
    fig_num=1;
end
% Plot state populations
for i=1:n_fluorophores
    my_fig=figure(fig_num);
    my_fig.Name=sm_par(i).fluorophore_name;
    clf
    hold on
    states_to_plot=unique(sm_par(i).state_ids);
    for k=1:numel(states_to_plot)
        plot(1:n_images,ens(i).det_p(states_to_plot(k),:),'LineWidth',3)
    end
    xlabel('Frame #');
    ylabel('Normalized Population');
    my_title=strrep(smis_title,'_','\_');
    my_fluoname=strrep(sm_par(i).fluorophore_name,'_','\_');
    title([my_fluoname, ' ensemble state population: ',my_title]);
    set(gca,'FontSize',16)
    set(gca,'LineWidth',2)
    my_legend=ens(i).state_names(states_to_plot);
    for j=1:numel(states_to_plot)
        my_legend{j}=strrep(my_legend{j},'_','\_');
    end
    legend(my_legend)
    hold off
    fig_num=fig_num+1;
end

% Plot the signal
if im_par.two_channel==1
    n_channels=2;
else
    n_channels=1;
end

for i=1:n_fluorophores
    my_fig=figure(fig_num);
    my_fig.Name=sm_par(i).fluorophore_name;
    clf
    subplot(1,n_channels,1)
    hold on
    for k=1:sm_par(i).n_fluorescent_states
        plot(1:n_images,ens(i).S_ch1(k,:),'LineWidth',3)
    end
    xlabel('Frame #');
    ylabel('Fluorescence [# photons per molecule]');
    my_title=strrep(smis_title,'_','\_');
    my_fluoname=strrep(sm_par(i).fluorophore_name,'_','\_');
    title([my_fluoname, ' fluorescence signal (Channel 1): ',my_title]);
    set(gca,'FontSize',16)
    set(gca,'LineWidth',2)
    my_legend=ens(i).state_names(sm_par(i).fluorescent_states);
    for j=1:sm_par(i).n_fluorescent_states
        my_legend{j}=strrep(my_legend{j},'_','\_');
    end
    legend(my_legend)
    hold off
    
    if im_par.two_channel==1
        subplot(1,n_channels,2)
        hold on
        
        for k=1:sm_par(i).n_fluorescent_states
            plot(1:n_images,ens(i).S_ch2(k,:),'LineWidth',3)
        end
        
        xlabel('Frame #');
        ylabel('Fluorescence [# photons per molecule]');
        my_title=strrep(smis_title,'_','\_');
        title([my_fluoname, ' fluorescence signal (Channel 2): ',my_title]);
        set(gca,'FontSize',16)
        set(gca,'LineWidth',2)
        legend(my_legend)
        hold off
    end
    fig_num=fig_num+1;
end
