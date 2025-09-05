function [output_image] = add_frontier_zones(input_image, p, q)
    % ADD_FRONTIER_ZONES Adds frontier zones between existing regions in an image
    % Input:
    %   input_image - n by m matrix with values 0, 1, 2 representing zones
    %   p - thickness of frontier between central (2) and intermediate (1) zones
    %   q - thickness of frontier between intermediate (1) and peripheral (0) zones
    % Output:
    %   output_image - image with added frontier zones (values 3 and 4)
    
    % Validate input
    if ~ismatrix(input_image)
        error('Input must be a 2D matrix');
    end
    if ~all(ismember(unique(input_image), [0 1 2]))
        error('Input image must contain only values 0, 1, and 2');
    end
    
    % Initialize output image
    output_image = input_image;
    
    % Create frontier between central (2) and intermediate (1) zones
    if p > 0
        % Find boundaries between zone 2 and zone 1
        central_boundary = (input_image == 2) & imdilate(input_image == 1, strel('disk', 1));
        
        % Create dilated boundary with thickness p
        frontier_zone = imdilate(central_boundary, strel('disk', p-1));
        
        % Set frontier pixels to value 3 (only where not already zone 1 or 2)
        output_image(frontier_zone & ~(input_image == 1 | input_image == 2)) = 3;
    end
    
    % Create frontier between intermediate (1) and peripheral (0) zones
    if q > 0
        % Find boundaries between zone 1 and zone 0
        intermediate_boundary = (input_image == 1) & imdilate(input_image == 0, strel('disk', 1));
        
        % Create dilated boundary with thickness q
        frontier_zone = imdilate(intermediate_boundary, strel('disk', q-1));
        
        % Set frontier pixels to value 4 (only where not already zone 0 or 1)
        output_image(frontier_zone & ~(input_image == 0 | input_image == 1)) = 4;
    end
end