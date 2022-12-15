function [simulation_files,continue_OK]=prepare_files(outfiledir,outfilename,lasers,force_overwrite)

n_lasers=size(lasers,2);

simulation_files=struct(...
    'im_ch1_out',fullfile(outfiledir,[outfilename,'_ch1.tif']), ...
    'im_ch2_out',fullfile(outfiledir,[outfilename,'_ch2.tif']), ...
    'data_out',fullfile(outfiledir,[outfilename,'.mat']), ...
    'diary_out',fullfile(outfiledir,[outfilename,'_diary.txt']), ...
    'im_true_out',fullfile(outfiledir,[outfilename,'_true.tif']), ...
    'im_true_out_3D',fullfile(outfiledir,[outfilename,'_true.mat']), ...
    'im_ch1_dl_out',fullfile(outfiledir,[outfilename,'_ch1_dl.tif']), ...
    'im_ch2_dl_out',fullfile(outfiledir,[outfilename,'_ch2_dl.tif']) ...
    );

% Add lasers filenames
for k=1:n_lasers
    simulation_files.(lasers(k).name)=fullfile(outfiledir,[outfilename,'_',lasers(k).name,'.tif']);
end

fields = fieldnames(simulation_files);

files_exist=0;
for k=1:size(fields,1)
    if exist(simulation_files.(fields{k}),'file')
        files_exist=1;
    end
end

if files_exist==1 && force_overwrite==0
    fig = uifigure;
    fig.Position = [800 600 400 160];
    MySel = uiconfirm(fig,'Simulation files already exist ! Overwrite ?','','Icon','warning');
    close(fig)
    
    diary off
    if strcmp(MySel,'OK')
        for k=1:size(fields,1)
            if exist(simulation_files.(fields{k}),'file')
                delete(simulation_files.(fields{k}));
            end
        end
        continue_OK=1;
    else
        continue_OK=0;
        return;
    end
elseif force_overwrite==1
    diary off
    for k=1:size(fields,1)
        if exist(simulation_files.(fields{k}),'file')
            delete(simulation_files.(fields{k}));
        end
    end
    continue_OK=1;
else
    continue_OK=1;
end
