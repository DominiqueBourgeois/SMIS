function [x,y] = wheresub(subscript,array)
% NAME:
%	WHERESUB
%
% PURPOSE:
%	Convert a 1d subscript value or array of value to a pixel position 
%	in the 2d array, or to an array containing the ensemble of these
%	pixel position
%
% CATEGORY:
%	Signal, image processing.
%
% CALLING SEQUENCE:
%	xy = wheresub(subscript,array)
%
% INPUTS:
%	subscript: the subscript or array of subscripts to process
%	array:	The array where to locate the subscript
%
% OUTPUTS:
%	x = the x coordinates
%   y = the y coordinates
%
% COMMON BLOCKS:
%	None.
%
% SIDE EFFECTS:
%	None.
%
% RESTRICTIONS:
%	None.
%
% MODIFICATION HISTORY:
%	D.Bourgeois, April 2011.

%get the size of the image array
size_im=size(array);

%indices in Matlab run column after column
 y = fix((subscript-1)/size_im(1))+1; % x is the row # in im
 x = subscript - (y-1)*size_im(1); % y is the column # in im
end

