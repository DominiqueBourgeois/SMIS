function im_par=init_EMCCD_frames(im_par,outfiledir,outfilename)
%
% PURPOSE:
%   Init Tiff Files containing EMCCD frames along acquisition
%
% INPUTS:
%	im_par: the imaging parameters
%   outfiledir: Directory of output file
%   outfilename: Filename four output file
%
% OUTPUTS:
%	im_par: the updated imaging parameters
%
% MODIFICATION HISTORY:
%	D.Bourgeois, December 2021.

%Estimate the size of the output stacks
StackSize=im_par.n*im_par.m*im_par.n_images*2e-9; % Image Size in GB
if StackSize>=4
    BigTif=1; % Set to 1 to write big tiff images > 4MB  (more tedious to open with image J)
else
    BigTif=0;
end

%Open the Tiff file at the first frame
if BigTif==0
    MyTiffFile_Ch1 = Fast_Tiff_Write(fullfile(outfiledir,[outfilename,'_ch1.tif']),1,0);
else
    MyTiffFile_Ch1 = Fast_BigTiff_Write(fullfile(outfiledir,[outfilename,'_ch1.tif']),1,0);
end

if im_par.two_channel==1 && im_par.single_CCD==0 % In that case output a ch2 image
    if BigTif==0
        MyTiffFile_Ch2 = Fast_Tiff_Write(fullfile(outfiledir,[outfilename,'_ch2.tif']),1,0);
    else
        MyTiffFile_Ch2 = Fast_BigTiff_Write(fullfile(outfiledir,[outfilename,'_ch2.tif']),1,0);
    end    
else
    MyTiffFile_Ch2=[];
end

im_par.TiffInfo.MyTiffFile_Ch1=MyTiffFile_Ch1;
im_par.TiffInfo.MyTiffFile_Ch2=MyTiffFile_Ch2;
