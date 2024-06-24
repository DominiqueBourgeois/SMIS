function overlay_SMIS_tracks_on_data(tracks, stack, im_par, par)

% PURPOSE:
% See how tracks match with raw data
%
% INPUTS:
%	tracks: the tracks
%   stack: the stack containing the raw data
%   im_par: the imaging parameters
%   par: additional parameters
%
% OUTPUTS:
%
% MODIFICATION HISTORY:
%	D.Bourgeois, January 2022


close all

raster=im_par.raster; % [nm]

track_range=par.track_range; %track range to analyse

plot_tracks=par.plot_tracks; % plot X% of the selected tracks

plot_tracks_on_image=1; % Interactively plot the tracks superimposed on frames
%if yes load the image stack
color_track_mode=1; % 1: One color per track; 0: Color gradient per track
overlay_full_track=1; % 1: Full track overlaid on image; 0: Only track until current frame is shown
% show_matched_tracks=2; % 0: Do not show, 1: Show only main matched track; 2: Show all matched track;
area_size=20; % minimum area size for display


%% END OF INPUT

s=size(stack);
%Converting tracks to localizations

GT_tracks=tracks;
GT_ntracks=length(GT_tracks);

GT_tl=[GT_tracks.tracklength];
GT_tracks=GT_tracks(GT_tl >= track_range(1) & GT_tl <= track_range(2));
GT_ntracks_selected=length(GT_tracks);

GT_locs=Ground_truth_tracks_to_localizations(GT_tracks,raster);


if plot_tracks>0
    figure(1);
    set(gcf, 'WindowStyle', 'docked')
    clf; set(gcf, 'color', [1 1 1]);
    %     set(gcf, 'Position', get(0, 'Screensize'));
    hold on;
    axis image
    for k=1:GT_ntracks_selected
        if k*plot_tracks/100==fix(k*plot_tracks/100)
            clc; disp(['Tracks [%] : ',num2str(100*k/GT_ntracks_selected)]);
            
            GT_t=GT_tracks(k); % GT tracks to be analyzed
            
            ColorVector = rand(1,3);
            plot([GT_t.x], [GT_t.y], '--+', 'Color', ColorVector, 'LineWidth',2);
        end
    end
end


if plot_tracks_on_image>0 && plot_tracks>0
    
    n_frames=par.n_frames;
    
    disp_fig=figure;
    
    set(disp_fig, 'WindowStyle', 'docked')
    colormap('gray');
    contrast_factor=1; % to modulate contrast
    
    X_all=[GT_locs.x]/raster; %
    Y_all=[GT_locs.y]/raster;
    
    
    ColorVector_db = rand(100,3);
    while 1
        figure(1)
        disp('Select spot on track image ...');
        [x,y] = ginput(1);
        %select closest track
        d2=(x-X_all).^2+(y-Y_all).^2;
        k=find(d2==min(d2),1);
        GT_t=GT_tracks(GT_locs(k).tn); % GT tracks to be analyzed
        
        t='n'; % String to select frame
        
        GT_x=[GT_t.x];
        GT_y=[GT_t.y];
        GT_f=[GT_t.frames];
        
        disp(['Frame range: from: ', num2str(min(GT_f)), ' to: ', num2str(max(GT_f))]);
        
        %Show the frame
        frame_to_show=GT_f(1);
        
        min_x=max(1,fix(min(GT_x)-area_size/2)); min_y=max(1,fix(min(GT_y)-area_size/2));
        max_x=min(s(2),ceil(max(GT_x)+area_size/2)); max_y=min(s(1),ceil(max(GT_y)+area_size/2));
        
        figure(disp_fig)
        while t~='s'
            if strcmp(t,'f')==1 ; frame_to_show=min([frame_to_show+1,n_frames]); end
            if strcmp(t,'b')==1 ; frame_to_show=max([frame_to_show-1,1]); end
            if strcmp(t,'ff')==1
                c_frame_to_show=frame_to_show;
                frame_to_show=str2double(input('Which frame to look at ? ','s'));
                if frame_to_show<1 || frame_to_show>n_frames
                    disp('Cannot reach this frame!');
                    frame_to_show=c_frame_to_show;
                end
            end
            
            area=stack(min_y:max_y,min_x:max_x, frame_to_show);
            area_r=imrotate(area,90); % to match track display on figure 1
            
            clims=1/contrast_factor*[0.8*min(area(:)) 1.2*max(area(:))];
            imagesc(area_r,clims); axis image;
            
            %show the track(s)
            hold on
            
            
            if overlay_full_track==1
                X=GT_x-min_x+1; Y=GT_y-min_y+1;
            else
                X=GT_x(GT_f<=frame_to_show)-min_x+1; Y=GT_y(GT_f<=frame_to_show)-min_y+1;
            end
            Xr=Y; Yr=size(area,2)+1-X; %rotate coordinates to match figure 1
            
            if color_track_mode==1
                ColorVector = ColorVector_db(end,:);
                plot(Xr, Yr, '--+', 'Color', ColorVector, 'LineWidth',2);
            elseif color_track_mode==0 %  color gradient per track
                ColorVector = ColorCoding(length(X));
                for n = 1 : length(X)-1
                    plot([Xr(n),Xr(n+1)], [Yr(n),Yr(n+1)], '--+', 'Color', ColorVector(n, :), 'LineWidth',2);
                end
                plot(Xr(end), Yr(end), '+', 'Color', ColorVector(end, :), 'LineWidth',2);
            end
            
            hold off
            
            t=input('Type ''f'' to move forward, ''b'' backward, ''ff'' to reach a frame, ''i'' to change contrast, ''n'' to next track, ''q'' to Quit: ','s');
            if t=='i'
                disp(['Actual contrast factor: ', num2str(contrast_factor)]);
                contrast_factor=str2double(input('Enter contrast factor: ','s'));
            end
            clc % Clear screen
            if t=='n'; break; end
            if t=='q'; return; end
        end
    end
end

% if reference_tracks==1
disp(['Number of GT tracks: ', num2str(GT_ntracks)]);
disp(['Number of GT tracks within selected range (', num2str(track_range), '): ' num2str(GT_ntracks_selected)]);
disp(['Number of SW tracks: ', num2str(SW_ntracks)]);

