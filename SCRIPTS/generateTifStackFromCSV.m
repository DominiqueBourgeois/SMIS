function generateTifStackFromCSV(filename, R, selected_ids)
    % Read the CSV file
    data = readtable(filename);

    % Extract X, Y, Z columns
    X = data.x;
    Y = data.y;
    Z = data.z;
    ids = data.particle+1; % Starts at 0, so add +1


    % Select the elements that match
    logicalIndex = ismember(ids, selected_ids);
    X=X(logicalIndex);
    Y=Y(logicalIndex);
    Z=Z(logicalIndex);
    
    % Convert coordinates to indices based on resolution R
    X_idx = round(X / R);
    Y_idx = round(Y / R);
    Z_idx = round(Z / R);

    % Make sure coodinates are positive
    X_idx=X_idx-min(X_idx)+1;
    Y_idx=Y_idx-min(Y_idx)+1;
    Z_idx=Z_idx-min(Z_idx)+1;

    % Determine the size of the 3D matrix
    maxX = max(X_idx);
    maxY = max(Y_idx);
    maxZ = max(Z_idx);

    % Initialize a 3D matrix of zeros
    tifStack = zeros(maxX, maxY, maxZ);

    % Set the voxels corresponding to the points to 1
    for i = 1:length(X_idx)
        tifStack(X_idx(i), Y_idx(i), Z_idx(i)) = 1;
    end

    % Write the 3D matrix to a .tif file
    for z = 1:maxZ
        if z == 1
            imwrite(tifStack(:, :, z), 'output_stack.tif', 'Compression', 'none');
        else
            imwrite(tifStack(:, :, z), 'output_stack.tif', 'WriteMode', 'append', 'Compression', 'none');
        end
    end

    % Display the 3D tif stack interactively
    displayTifStack('output_stack.tif', maxZ);
end

function displayTifStack(filename, maxZ)
    % Read the .tif stack
    tifStack = zeros(size(imread(filename)));
    for z = 1:maxZ
        tifStack(:, :, z) = imread(filename, z);
    end

    % Display the stack interactively
    figure;
    vol3d('cdata',tifStack);
    colormap('jet');
    axis image
    view(3);
    % for z = 1:maxZ
    %     imagesc(tifStack(:, :, z));
    %     title(['Slice ', num2str(z)]);
    %     colormap(gray);
    %     axis image;
    %     drawnow;
    %     pause(0.1);
    % end
end
