function [subimages, n_frames] = read_SMIS_MultipageTiff(filename, frame_id)

% Read a multipage tiff, assuming each file is the same size
n_frames= size(imfinfo(filename),1);
t = Tiff(filename, 'r');

if nargin > 1 % read specific frame
    setDirectory(t, frame_id);
    subimages = t.read();
   
else
    
    subimages(:,:,1) = t.read(); % Read the first image to get the array dimensions correct.
    if t.lastDirectory()
        return; % If the file only contains one page, we do not need to continue.
    end
    
    % Read all remaining pages (directories) in the file
    t.nextDirectory();
    i=0;
    warning('off','imageio:tiffmexutils:libtiffWarning')
    while true
        
        if ~mod(i,100)
            clc;
            disp(['Loading stack ... ',num2str(100*i/n_frames), '%']);
        end
        subimages(:,:,end+1) = t.read();
        if t.lastDirectory()
            break;
        else
            t.nextDirectory();
            i=i+1;
        end
    end
     disp('Done !');
end
end


