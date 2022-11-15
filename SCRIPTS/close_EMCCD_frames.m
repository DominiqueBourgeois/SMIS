function close_EMCCD_frames(im_par)
%
% PURPOSE:
%   Close Tiff Files containing EMCCD frames along acquisition
%
% INPUTS:
%	im_par: the imaging parameters
%
% OUTPUTS:

% MODIFICATION HISTORY:
%	D.Bourgeois, December 2021.

im_par.TiffInfo.MyTiffFile_Ch1.close;
if im_par.two_channel==1 && im_par.single_CCD==0 % In that case output a ch2 image
    im_par.TiffInfo.MyTiffFile_Ch2.close;
end
