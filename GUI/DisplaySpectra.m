function spectra_figure_number=DisplaySpectra(smis_gui_parameters)
%
% PURPOSE:
%   Display fluorophore spectra
%
% INPUTS:
% smis_gui_parameters: the GUI parameters
%
% OUTPUTS:
%   spectra_figure_number: updated spectra_figure_number
%
% MODIFICATION HISTORY:
%	D.Bourgeois, February 2021.
%	D.Bourgeois, March 2021. Correct for proper handling of figure creation


% Create a figure if necessary
create_figure=1;
% The following code ensures that the right figure is addressed
% Indeed, it might be that the figure has been created before, then deleted, and the same figure number has been re-created for another purpose and is still open

% Take a look at the current figures open
figHandles = findobj('Type', 'figure'); % Get the array of all figures open
n_fig=size(figHandles,1);
for k=1:n_fig
    if strcmp(figHandles(k).Name,'Abs/Exc and Emission Spectra')==1 % Our figure already exists
        create_figure=0;
        spectra_figure_number=figHandles(k).Number;
    end
end



% create a figure if necessary
if create_figure==1
    scrsz=get(0,'ScreenSize');
    w_xsize=scrsz(3)*3/4;
    w_ysize=scrsz(4)*3/4;
    spectra_figure=figure;
    spectra_figure_number=spectra_figure.Number;

    colormap(gray);
    set(spectra_figure,'color','white');
    set(spectra_figure,'Position',[10 scrsz(4)/10 w_xsize w_ysize])
    set(spectra_figure,'color','white');
    set(spectra_figure,'Name','Abs/Exc and Emission Spectra');
    set(spectra_figure,'NumberTitle','off');
else
    figure(spectra_figure_number);
    clf
end


% Extract all relevant parameters
Fluorophores=smis_gui_parameters.Fluorophores;
TwoChannel=smis_gui_parameters.TwoChannel;

% Get the number of fluorophores
n_fluorophores=size(Fluorophores,2);

%Figure out the number of columns
n_spectra=zeros(1,n_fluorophores);
for k=1:n_fluorophores
    n_spectra(k)=numel(Fluorophores(k).Photophysics.photoactive_states);
end

% Display differently if only one dye or several
if n_fluorophores==1
    n_rows=2;
    n_columns=Fluorophores(1).Photophysics.n_fluorescent_states; % Display the Fluorescent states on the top
else
    n_rows=n_fluorophores; % Display all spectra for one dye on a single row
    n_columns=max(n_spectra);
end

% MyFontsize=fix(18*4/(n_rows*n_columns));
% MyFontsize=max([MyFontsize,7]);
% MyFontsize=min([MyFontsize,30]);

MyFontsize=fix(5*4/(n_rows*n_columns));
MyFontsize=max([MyFontsize,7]);
MyFontsize=min([MyFontsize,20]);

MyLineWidth=fix(1*4/(n_rows*n_columns));
MyLineWidth=max([MyLineWidth,2]);
MyLineWidth=min([MyLineWidth,4]);


% Go through the fluorophores
for k=1:n_fluorophores
    Photophysics=Fluorophores(k).Photophysics;
    % First plot spectra for the fluorescence states: absorbance/excitation and
    % emission
    for j=1:Photophysics.n_fluorescent_states

        state_name=Photophysics.state_names{Photophysics.fluorescent_states(j)};
        subplot(n_rows, n_columns, (k-1)*n_columns+j);
        hold on
        %First Plot excitation spectrum
        exc_spectrum=Photophysics.exc_spectra(j).s;
        tmp_plot=plot(exc_spectrum(:,1),exc_spectrum(:,2));
        % Get the right color for the spectrum
        % get the Peak maximum Above 300 nm (center of mass does not work well, prevent Selecting 260 nanometer peak)
        sub_s=exc_spectrum(exc_spectrum(:,1)>300,:);
        pm=sub_s(sub_s(:,2)==max(sub_s(:,2)),1);
        sc=wavelength2rgb(-pm(1)); % Negative value has to be used for a single-point entry
        set(tmp_plot,'LineWidth',1.5);
        set(tmp_plot,'Color',sc(1,:));

        %Second, Plot emission spectrum
        em_spectrum=Photophysics.em_spectra(j).s;
        tmp_plot=plot(em_spectrum(:,1),em_spectrum(:,2));
        % Get the right color for the spectrum
        % get the Peak maximum Above 300 nm (center of mass does not work well, prevent Selecting 260 nanometer peak)
        sub_s=em_spectrum(em_spectrum(:,1)>300,:);
        pm=sub_s(sub_s(:,2)==max(sub_s(:,2)),1);
        sc=wavelength2rgb(-pm(1)); % Negative value has to be used for a single-point entry
        set(tmp_plot,'LineWidth',1.5);
        set(tmp_plot,'Color',sc(1,:));

        ylim([0 max([max(em_spectrum(:,2)),max(exc_spectrum(:,2))])]);
        xlim([min(exc_spectrum(:,1)) max(em_spectrum(:,1))]);
        xlabel('Wavelength [nm]','fontsize',MyFontsize,'fontweight','b')
        ylabel('Exc/Em [AU]','fontsize',MyFontsize,'fontweight','b')

        %         MyTitle=check_underscore([Fluorophores(k).name ': Exc/Em spectrum of state: ', newline, state_name]);
        MyTitle=check_underscore([Fluorophores(k).name, newline, state_name]);
        title(gca,MyTitle,'FontSize',MyFontsize,'FontWeight','bold');

        % now Plot the filters
        if TwoChannel.use_exp_ch1_filter==0
            [x_ch1,y_ch1,~] = get_filter_profile(eval(TwoChannel.em_filter_1));
        else
            x_ch1=TwoChannel.exp_ch1_filter(1,:);
            y_ch1=TwoChannel.exp_ch1_filter(2,:);
        end
        tmp_plot=plot(x_ch1,y_ch1,'--'); %Plot filter1 in dashed line
        set(tmp_plot,'LineWidth',1.5);
        set(tmp_plot,'Color','green');
        if TwoChannel.state==1
            if TwoChannel.use_exp_ch2_filter==0
                [x_ch2,y_ch2,~] = get_filter_profile(eval(TwoChannel.em_filter_2));
            else
                x_ch2=TwoChannel.exp_ch2_filter(1,:);
                y_ch2=TwoChannel.exp_ch2_filter(2,:);
            end
            tmp_plot=plot(x_ch2,y_ch2,'--'); %Plot filter1 in dashed line
            set(tmp_plot,'LineWidth',1.5);
            set(tmp_plot,'Color','red');

            % Eventually plot dichroic
            if TwoChannel.add_dichroic==1
                if TwoChannel.use_exp_dichroic==0
                    [x_df,y_df,~] = get_dichroic_filter_profile(eval(TwoChannel.dichroic_filter));
                else
                    x_df=TwoChannel.exp_dichroic_filter(1,:);
                    y_df=TwoChannel.exp_dichroic_filter(2,:);
                end
                tmp_plot=plot(x_df,y_df,'-.'); %Plot dichroic filter in dashed dotted line
                set(tmp_plot,'LineWidth',1.5);
                set(tmp_plot,'Color','blue');
            end

        end

        set(gca,'LineWidth',MyLineWidth)
        set(gca,'FontSize',MyFontsize);
        set(gca,'FontWeight','b');


        hold off;
    end

    % Then plot Absorbance spectra for the Photoactive dark states states
    for j=1:Photophysics.n_photoactive_dark_states
        photoactive_dark_states=setdiff(Photophysics.photoactive_states,Photophysics.fluorescent_states);
        state_name=Photophysics.state_names{photoactive_dark_states(j)};

        if n_fluorophores==1
            subplot(n_rows, Photophysics.n_photoactive_dark_states, j+Photophysics.n_photoactive_dark_states);
        else
            subplot(n_rows, n_columns, (k-1)*n_columns+j+Photophysics.n_fluorescent_states);
        end
        hold on
        % Plot Absorbance spectrum
        abs_spectrum=Photophysics.dark_spectra(j).s;
        tmp_plot=plot(abs_spectrum(:,1),abs_spectrum(:,2));
        % Get the right color for the spectrum
        % get the Peak maximum Above 300 nm (center of mass does not work well, prevent Selecting 260 nanometer peak)
        sub_s=abs_spectrum(abs_spectrum(:,1)>300,:);
        pm=sub_s(sub_s(:,2)==max(sub_s(:,2)),1);
        sc=wavelength2rgb(-pm(1)); % Negative value has to be used for a single-point entry
        set(tmp_plot,'LineWidth',1.5);
        set(tmp_plot,'Color',sc(1,:));

        ylim([0 max(abs_spectrum(:,2))]);
        xlim([min(abs_spectrum(:,1)) max(abs_spectrum(:,1))]);
        xlabel('Wavelength [nm]','fontsize',MyFontsize,'fontweight','b')
        ylabel('Abs [AU]','fontsize',MyFontsize,'fontweight','b')
        %         MyTitle=check_underscore([Fluorophores(k).name ': Dark spectrum of state: ', newline, state_name]);
        MyTitle=check_underscore([Fluorophores(k).name, newline, state_name]);
        title(gca,MyTitle,'FontSize',MyFontsize,'FontWeight','bold');

        set(gca,'LineWidth',MyLineWidth)
        set(gca,'FontSize',MyFontsize);
        set(gca,'FontWeight','b');

        hold off;
    end
end

