function psf_figure_number=DisplayPsf(smis_gui_parameters)
%
% PURPOSE:
%   Display psfs shape
%
% INPUTS:
%	smis_gui_parameters: the gui parameters
%
% OUTPUTS:
%   psf_figure_number: updated psf_figure_number
%
% MODIFICATION HISTORY:
%	D.Bourgeois, February 2021.
%	D.Bourgeois, March 2021. Correct for proper handling of figure creation

% Extract all relevant parameters
Fluorophores=smis_gui_parameters.Fluorophores;
TwoChannel=smis_gui_parameters.TwoChannel;
obj_na=smis_gui_parameters.obj_and_psf.obj_na;
psf_figure_number=smis_gui_parameters.psf_figure_number;

simul_3D=smis_gui_parameters.simul_3D;

if simul_3D==1 && isempty(smis_gui_parameters.image_size.nz)
    warndlg('3D Virtual sample must be loaded first !','Warning')
    return
end

% Get the number of fluorophores
n_fluorophores=size(Fluorophores,2);

% First calculate the PSF width
% First go with channel 1

%Produce temporary im_par and sm_par to use the same routines as in
%smis_main
im_par.filters=smis_gui_parameters.TwoChannel;
im_par.filters.ch1_filter=eval(im_par.filters.em_filter_1);
im_par.filters.ch2_filter=eval(im_par.filters.em_filter_2);
im_par.two_channel=smis_gui_parameters.TwoChannel.state;
im_par.obj.na=smis_gui_parameters.obj_and_psf.obj_na;
im_par.raster=smis_gui_parameters.raster;
im_par.two_channel_defocus=smis_gui_parameters.TwoChannel.defocus;

if simul_3D==1
    im_par.psf_n_zslices=smis_gui_parameters.obj_and_psf.psf_n_zslices;
    im_par.depth_of_focus=smis_gui_parameters.obj_and_psf.obj_depth_of_focus;
    im_par.sample_zcenter=smis_gui_parameters.obj_and_psf.sample_zcenter;
    im_par.nz=smis_gui_parameters.image_size.nz;
    im_par.psf_astigmatism_x=smis_gui_parameters.obj_and_psf.psf_astigmatism_x;
    im_par.psf_astigmatism_y=smis_gui_parameters.obj_and_psf.psf_astigmatism_y;
end

sm_par(1:n_fluorophores)=struct;
for i=1:n_fluorophores
    sm_par(i).n_fluorescent_states=Fluorophores(i).Photophysics.n_fluorescent_states;
    sm_par(i).spectral_data.em_spectra=Fluorophores(i).Photophysics.em_spectra;
    sm_par(i).spectral_data.n_fluorescent_states=sm_par(i).n_fluorescent_states;
    sm_par(i).fluorophore_name=Fluorophores(i).name;
    sm_par(i).filter_profiles=[];
end

%get filter profiles
for i=1:n_fluorophores
    sm_par(i) = get_sm_filter_profiles(im_par,sm_par(i));
end

%get PSF width
sm_par=get_psf_width(n_fluorophores,sm_par,im_par);

% get PSFs
psf_mode='Simple_Gaussian';
sm_par=get_psf(n_fluorophores, sm_par, im_par, psf_mode, 1, simul_3D);
if im_par.two_channel==1 % Set psf for 2nd channel
    sm_par=get_psf(n_fluorophores, sm_par, im_par, psf_mode,2, simul_3D);
end

% Create a figure if necessary
create_figure=1;
% The following code ensures that the right figure is addressed
% Indeed, it might be that the figure has been created before, then deleted, and the same figure number has been re-created for another purpose and is still open

% Take a look at the current figures open
figHandles = findobj('Type', 'figure'); % Get the array of all figures open
n_fig=size(figHandles,1);
for k=1:n_fig
    if strcmp(figHandles(k).Name,'Point Spread Functions')==1 % Our figure already exists
        create_figure=0;
        psf_figure_number=figHandles(k).Number;
    end
end


% create a figure if necessary
if create_figure==1
    scrsz=get(0,'ScreenSize');
    w_xsize=scrsz(3)*3/4;
    w_ysize=scrsz(4)*3/4;
    psf_figure=figure;
    psf_figure_number=psf_figure.Number;

    colormap(gray);
    set(psf_figure,'color','white');
    set(psf_figure,'Position',[10 scrsz(4)/10 w_xsize w_ysize])
    set(psf_figure,'color','white');
    set(psf_figure,'Name','Point Spread Functions');
    set(psf_figure,'NumberTitle','off');
else
    figure(psf_figure_number);
    clf;
end

n_fluorescent_states=zeros(1,n_fluorophores);
for i=1:n_fluorophores
    n_fluorescent_states(i)=sm_par(i).n_fluorescent_states;

end

n_rows=sum(n_fluorescent_states); % Number of rows
if im_par.two_channel==1
    n_columns=2; % Number of columns
else
    n_columns=1; % Number of columns
end

MyFontsize=min([fix(18*4/(n_rows*n_columns)),20]);
MyFontsize=max([MyFontsize,7]);
MyFontsize=min([MyFontsize,30]);

% Display the psf profiles in channel 1
row_n=1;
for k=1:n_fluorophores
    fluorescent_states=Fluorophores(k).Photophysics.fluorescent_states;
    for j=1:n_fluorescent_states(k)
        psf=sm_par(k).psf_par_ch1(j).psf;
        if simul_3D==0
            if ~any(isnan(psf)) % Check that the PSF does exist
                subplot(n_rows,n_columns,row_n)
                imagesc(psf, [0,max(max(psf))]);
                axis image;
                xlabel('X Coordinate [raster]','fontsize',MyFontsize,'fontweight','b')
                ylabel('Y Coordinate [raster]','fontsize',MyFontsize,'fontweight','b')
            end
        elseif simul_3D==1
            if ~any(isnan(psf(1).plane)) % Check that the PSF does exist
                subplot(n_rows,n_columns,row_n)
                ShowPSF3D(psf,MyFontsize);
            end
        end
        state_name=Fluorophores(k).Photophysics.state_names{fluorescent_states(j)};

        title(gca,['Channel 1 : ',strrep(Fluorophores(k).name,'_',' '),', ', state_name], 'FontWeight','bold','FontSize',MyFontsize);
        row_n=row_n+1;

        % Display the psf profiles in channel 1
        if im_par.two_channel==1
            psf=sm_par(k).psf_par_ch2(j).psf;
            if simul_3D==0
                if ~any(isnan(psf)) % Check that the PSF does exist
                    subplot(n_rows,n_columns,row_n)
                    imagesc(psf, [0,max(max(psf))]);
                    axis image;
                    xlabel('X Coordinate [raster]','fontsize',MyFontsize,'fontweight','b')
                    ylabel('Y Coordinate [raster]','fontsize',MyFontsize,'fontweight','b')
                end
            elseif simul_3D==1
                if ~any(isnan(psf(1).plane)) % Check that the PSF does exist
                    subplot(n_rows,n_columns,row_n)
                    ShowPSF3D(psf,MyFontsize);
                end
            end
            title(gca,['Channel 2 : ',strrep(Fluorophores(k).name,'_',' '),', ', state_name], 'FontWeight','bold','FontSize',MyFontsize);
            row_n=row_n+1;
        end
    end
end
end



function ShowPSF3D(psf,MyFontsize)
%Limit to 32 pixels
size_lim=32;
rbox_x_max=max([psf.rbox_x]);
rbox_y_max=max([psf.rbox_y]);
rbox_z=size(psf,2);
psf3D=zeros(2*rbox_y_max+1,2*rbox_x_max+1,rbox_z);

for i=1:rbox_z
    psf3D(rbox_y_max-psf(i).rbox_y+1:rbox_y_max+psf(i).rbox_y+1,rbox_x_max-psf(i).rbox_x+1:rbox_x_max+psf(i).rbox_x+1,i)=psf(i).plane;
end
if size(psf3D,1)>size_lim
    cx=fix(size(psf3D,1)/2.0);
    sl=fix(size_lim/2.0);
    psf3D([1:cx-sl,cx+sl+1:end],:,:)=[];
end
if size(psf3D,2)>size_lim
    cy=fix(size(psf3D,2)/2.0);
    sl=fix(size_lim/2.0);
    psf3D(:,[1:cy-sl,cy+sl+1:end],:)=[];
end

vol3d('cdata', psf3D);
colormap(gray);
axis image
xlabel('X Coordinate [raster]','fontsize',MyFontsize,'fontweight','b')
ylabel('Y Coordinate [raster]','fontsize',MyFontsize,'fontweight','b')
zlabel('Z Coordinate [raster]','fontsize',MyFontsize,'fontweight','b')
view(3);
end
