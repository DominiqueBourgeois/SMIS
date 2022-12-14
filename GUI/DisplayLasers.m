
function lasers_figure_number=DisplayLasers(smis_par)
%
% PURPOSE:
%   Display lasers shape and intensity profile from the SMIS GUI
%
% INPUTS:
%	lasers: the lasers
%	n_images: the number of images
%   lasers_figure_number: Handle to the figure number displaying lasers
%
% OUTPUTS:
%   lasers_figure_number: updated lasers_figure_number
%
% MODIFICATION HISTORY:
%	D.Bourgeois, February 2021.
%	D.Bourgeois, March 2021. Correct for proper handling of figure creation
%	D.Bourgeois, September 2022. Add raster and tidy up figure


lasers=smis_par.lasers;
n_images=smis_par.n_images;
lasers_figure_number=smis_par.lasers_figure_number;
raster=smis_par.raster;

% Create a figure if necessary
create_figure=1;
% The following code ensures that the right figure is addressed
% Indeed, it might be that the figure has been created before, then deleted, and the same figure number has been re-created for another purpose and is still open

% Take a look at the current figures open
figHandles = findobj('Type', 'figure'); % Get the array of all figures open
n_fig=size(figHandles,1);
for k=1:n_fig
    if strcmp(figHandles(k).Name,'Lasers profiles and sequences')==1 % Our figure already exists
        create_figure=0;
        lasers_figure_number=figHandles(k).Number;
    end
end
    
if create_figure==1
    scrsz=get(0,'ScreenSize');
    w_xsize=scrsz(3)*3/4;
    w_ysize=scrsz(4)*3/4;
    lasers_figure=figure;
    lasers_figure_number=lasers_figure.Number;
    colormap(gray);
    set(lasers_figure,'color','white');
    set(lasers_figure,'Position',[10 scrsz(4)/10 w_xsize w_ysize])
    set(lasers_figure,'color','white');
    set(lasers_figure,'Name','Lasers profiles and sequences');
    set(lasers_figure,'NumberTitle','off');
else
    figure(lasers_figure_number);
    clf;
end

w_laser_ok=find([lasers.power]>0);
n_lasers=numel(w_laser_ok);


n_rows=n_lasers; % Number of rows
n_columns=2; % Number of columns

MyFontsize=fix(18*4/(n_rows*n_columns));
MyFontsize=max([MyFontsize,7]);
MyFontsize=min([MyFontsize,20]);

MyLineWidth=fix(3*4/(n_rows*n_columns));
MyLineWidth=max([MyLineWidth,2]);
MyLineWidth=min([MyLineWidth,5]);


% First display the laser profiles
for k=1:n_lasers
    if ~isempty(lasers(k).beam_profile) && lasers(k).power>0
        subplot(n_rows,n_columns,k)
        %Display laser profile with local laser power density 
        imagesc(lasers(k).beam_profile/((raster*1e-7)^2));
        axis image;
        title(gca,[check_underscore(lasers(k).name),' (',num2str(lasers(k).power),' mW)'], 'FontWeight','bold','FontSize',MyFontsize);
        xlabel('X [raster]')
        ylabel('Y [raster]')
        set(gca,'LineWidth',1);
        set(gca,'FontSize',MyFontsize);
        set(gca,'FontWeight','b');
    end
end



% Second display the laser sequences
for k=1:n_lasers
    if max(lasers(k).sequence*lasers(k).power)>0 && n_images>1
        if numel(lasers(k).wavelength)==1
            loc_c=wavelength2rgb(-lasers(k).wavelength); % Use negative value if numel(laser_wavelength)=1
        else
            loc_c=wavelength2rgb(lasers(k).wavelength);
        end
        
        %Display power density
        beam_power_density=lasers(k).sequence/100*lasers(k).max_beam_profile/((raster*1e-7)^2);
        subplot(n_rows,n_columns,n_rows+k)
%         scatter(1:n_images,lasers(k).sequence/100*lasers(k).power,[],loc_c,'filled')
        scatter(1:n_images,beam_power_density,[],loc_c,'filled')
        title(gca,[check_underscore(lasers(k).name),' laser sequence'], 'FontWeight','bold','FontSize',MyFontsize);
%         ylim([0 max(lasers(k).sequence/100*lasers(k).power)*1.2]);
        ylim([0 max(beam_power_density)*1.2]);
        xlabel('Frame #')
        ylabel('Power Density [W/cm??]')
        set(gca,'LineWidth',MyLineWidth)
        set(gca,'FontSize',MyFontsize);
        set(gca,'FontWeight','b');
    end
end
