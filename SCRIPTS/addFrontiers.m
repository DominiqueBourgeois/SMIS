function transformedImage = addFrontiers(image, p, q)
    % Create a copy of the original image
    transformedImage = double(image);

    % Create binary masks for each zone
    centralZone = (image == 2);
    intermediateZone = (image == 1);
    peripheralZone = (image == 0);

    % Define structuring elements for dilation
    seCentral = strel('disk', p);
    seIntermediate = strel('disk', q);

    % Dilate the central zone and subtract the original to get the frontier
    dilatedCentral = imdilate(centralZone, seCentral);
    frontierCentralIntermediate = dilatedCentral & intermediateZone;

    % Dilate the intermediate zone and subtract the original to get the frontier
    dilatedIntermediate = imdilate(intermediateZone, seIntermediate);
    frontierIntermediatePeripheral = dilatedIntermediate & peripheralZone;

    % Assign values to the frontiers in the original image
    transformedImage(frontierCentralIntermediate) = 3;
    transformedImage(frontierIntermediatePeripheral) = 4;

    % Ensure the output is an image matrix of type uint8
    transformedImage = uint8(transformedImage);
end

% Example usage:
% Define your image matrix 'image' here
% image = randi([0, 2], 100, 100); % Example random image
% p = 2; % Thickness of the frontier between central and intermediate zones
% q = 3; % Thickness of the frontier between intermediate and peripheral zones
% transformedImage = addFrontiers(image, p, q);
% imshow(transformedImage, []);
