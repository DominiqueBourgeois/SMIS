function [a_out,n,m,nz]=check_im_size(a_in, binning)
% NAME:
%	CHECK_IM_SIZE
%
% PURPOSE:
%	Crop image if binning wrong
%
% CATEGORY:
%
% CALLING SEQUENCE:
%
% INPUTS:
%
%
% OUTPUTS:
%

% MODIFICATION HISTORY:
%	D.Bourgeois, December 2011.

n_h=size(a_in,1);
m_h=size(a_in,2);
if numel(size(a_in))==3 % If 3D kernel
    nz_h=size(a_in,3);
    if (mod(n_h, binning)~=0) || (mod(m_h, binning)~=0) || (mod(nz_h, binning)~=0)
        crop_n=mod(n_h, binning);
        crop_m=mod(m_h, binning);
        crop_nz=mod(nz_h, binning);
        n_h=n_h-crop_n;
        m_h=m_h-crop_m;
        nz_h=nz_h-crop_nz;
        a_out=a_in(1:n_h,1:m_h,1:nz_h);
        disp(['Size of image (',num2str(size(a_in)),') divided by binning factor ',num2str(binning),' is non-integer: cropping image to (',num2str(size(a_out)),') !']);
%         warndlg(['Cropping image to (',num2str(size(a_out)),') !'],'Warning')
    else
        a_out=a_in;
    end
    n=fix(n_h/binning);
    m=fix(m_h/binning);
    nz=fix(nz_h/binning);
    
else % 2D mode
    if (mod(n_h, binning)~=0) || (mod(m_h, binning)~=0)
        crop_n=mod(n_h, binning);
        crop_m=mod(m_h, binning);
        n_h=n_h-crop_n;
        m_h=m_h-crop_m;
        a_out=a_in(1:n_h,1:m_h);
        disp(['Size of image (',num2str(size(a_in)),') divided by binning factor ',num2str(binning),' is non-integer: cropping image to (',num2str(size(a_out)),') !']);
%         warndlg(['Cropping image to (',num2str(size(a_out)),') !'],'Warning')
    else
        a_out=a_in;
    end
    n=fix(n_h/binning);
    m=fix(m_h/binning);
    nz=1;
end
end