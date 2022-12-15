function write_EMCCD_frames(im_par)
%
% PURPOSE:
%   Write EMCCD frames along acquisition: uses the Fast_Tiff writer 
%
% INPUTS:
%	im_par: the imaging parameters
%
% OUTPUTS:
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.
%	D.Bourgeois, December 2021.


% Save channel 1 and possibly channel 2 images
if max(max(im_par.emccd_im_ch1)) > (2^16-1)
    disp('WARNING! Channel 1 image contains saturated pixels: check parameters or reduce laser power !');
end

if im_par.two_channel==1 && im_par.single_CCD==1 % In that case output a single image that contains the two frames
    im_par.emccd_im_ch12(:,1:im_par.m)=im_par.emccd_im_ch1;
    im_par.emccd_im_ch12(:,im_par.m+1:end)=im_par.emccd_im_ch2;
    
    im_par.TiffInfo.MyTiffFile_Ch1.WriteIMG(uint16(im_par.emccd_im_ch12));
    
else % In that case output a ch1 image
    im_par.TiffInfo.MyTiffFile_Ch1.WriteIMG(uint16(im_par.emccd_im_ch1));
end
% writeDirectory(im_par.TiffInfo.MyTiffFile_Ch1);


if im_par.two_channel==1 && im_par.single_CCD==0 % In that case output a ch2 image
    if max(max(im_par.emccd_im_ch2)) > (2^16-1)
        disp('WARNING! Channel 2 image contains saturated pixels: check parameters or reduce laser power !');
    end
    im_par.TiffInfo.MyTiffFile_Ch2.WriteIMG(uint16(im_par.emccd_im_ch2));
end


