function write_EMCCD_frames_pct(det_im, im_par)
%
% PURPOSE:
%   Write EMCCD frames along acquisition: uses the Fast_Tiff writer 
%
% INPUTS:
%	det_im: the detector images 
%	im_par: the imaging parameters
%
% OUTPUTS:
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.
%	D.Bourgeois, December 2021.
%	D.Bourgeois, September 2022, introduce det_im.


% Save channel 1 and possibly channel 2 images
if max(max(det_im.emccd_im_ch1)) > (2^16-1)
    disp('WARNING! Channel 1 image contains saturated pixels: check parameters or reduce laser power !');
end

if im_par.two_channel==1 && im_par.single_CCD==1 % In that case output a single image that contains the two frames
    det_im.emccd_im_ch12(:,1:im_par.m)=det_im.emccd_im_ch1;
    det_im.emccd_im_ch12(:,im_par.m+1:end)=det_im.emccd_im_ch2;
    
    im_par.TiffInfo.MyTiffFile_Ch1.WriteIMG(uint16(det_im.emccd_im_ch12));
    
else % In that case output a ch1 image
    im_par.TiffInfo.MyTiffFile_Ch1.WriteIMG(uint16(det_im.emccd_im_ch1));
end
% writeDirectory(im_par.TiffInfo.MyTiffFile_Ch1);


if im_par.two_channel==1 && im_par.single_CCD==0 % In that case output a ch2 image
    if max(max(det_im.emccd_im_ch2)) > (2^16-1)
        disp('WARNING! Channel 2 image contains saturated pixels: check parameters or reduce laser power !');
    end
    im_par.TiffInfo.MyTiffFile_Ch2.WriteIMG(uint16(det_im.emccd_im_ch2));
end


