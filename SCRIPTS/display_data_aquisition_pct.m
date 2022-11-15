function display_par=display_data_aquisition_pct(det_im, im_par,display_par)
%
% PURPOSE:
%   Display EMCCD images during acquisition
%
% INPUTS:
%	det_im: the detector images 
%	im_par: the imaging parameters
%   display_par: The displaying parameters
%
% OUTPUTS:
%   display_par: The displaying parameters updated for figure number
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.
%	D.Bourgeois, September 2022, introduce det_im.

% create a figure
if (display_par.keep_current_figure==0 && im_par.current_frame==1) || isfield(display_par,'emccd_figure_number')==0 || ishandle(display_par.emccd_figure_number)==0
    w_xsize=display_par.scrsz(3)*3/4;
    w_ysize=display_par.scrsz(4)*1/2;
    display_par.emccd_figure_number=figure;colormap(gray);
    set(display_par.emccd_figure_number,'Position',[10 display_par.scrsz(4)/10 w_xsize w_ysize])
    set(display_par.emccd_figure_number,'color','white');
end
figure(display_par.emccd_figure_number);
fontsize=12;

n_rows=1; % Number of rows
if im_par.two_channel==1 % Number of columns
    n_columns=3;
else
    n_columns=2;
end

%First display TrueImage
%display "true" image as RGB
if im_par.current_frame==1
    subplot(n_rows,n_columns,1);
    if im_par.simul_3D==0
        imagesc(det_im.true_im_rgb);
        axis image;
    elseif im_par.simul_3D==1
        vol3d('cdata', det_im.true_3D_kernel);
        colormap('jet');
        axis image
        % axis equal off
        % set(gcf, 'color', 'w');
        view(3);
    end
    title(gca,'True image','FontWeight','bold','FontSize',fontsize);
    xlabel('X Coordinate [raster]','fontsize',fontsize,'fontweight','b')
    ylabel('Y Coordinate [raster]','fontsize',fontsize,'fontweight','b')
end


% Now display the EMCCD images
% First channel 1
if display_par.auto_display_range == 1
    display_par.disp_lim_ch1=[0,max(max(det_im.emccd_im_ch1))];
end
subplot(n_rows,n_columns,2);
imagesc(det_im.emccd_im_ch1, display_par.disp_lim_ch1);
axis image; colormap(gca,'gray')

title(gca,['Channel 1, Frame: ', num2str(im_par.current_frame)],...
    'FontWeight','bold','Color', [0 0 0]);
xlabel('X Coordinate [raster]','fontsize',fontsize,'fontweight','b')
ylabel('Y Coordinate [raster]','fontsize',fontsize,'fontweight','b')

% Second channel 2
if im_par.two_channel==1
    if display_par.auto_display_range == 1
        display_par.disp_lim_ch2=[0,max(max(det_im.emccd_im_ch2))];
    end
    subplot(n_rows,n_columns,3);
    imagesc(det_im.emccd_im_ch2, display_par.disp_lim_ch2);
    axis image; colormap(gca,'gray')
    
    title(gca,['Channel 2, Frame: ', num2str(im_par.current_frame)],...
        'FontWeight','bold','Color', [0 0 0]);
    xlabel('X Coordinate [raster]','fontsize',fontsize,'fontweight','b')
    ylabel('Y Coordinate [raster]','fontsize',fontsize,'fontweight','b')
end
pause(0.01);
