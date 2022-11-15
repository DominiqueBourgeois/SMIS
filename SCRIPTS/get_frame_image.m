function [sms,im_par]=get_frame_image(sms, sm_par, im_par)

%
% PURPOSE:
%   Get images as recorded by the detector
%
% INPUTS:
%   sms: the single molecules (with coordinates on high-resolution image) in raster units
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   im_par: updated for emccd_im_ch1 (image on channel 1) and emccd_im_ch2
%   (image on channel 2)
%   sms: The updated sms
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.
%	D.Bourgeois, May 2021. Introduce sm.n_phot_det_ch1 and sm.n_phot_det_ch2 in update_image and update_sub_image;

for k=1:sm_par.n_mol_eff
    if sms.sm(k).bleached==0 || sms.sm(k).bleached==im_par.current_frame
        if im_par.add_diffusion==0 || im_par.use_diffuse_psf==0 || isempty(sms.sm(k).sub_x) % general case
            [x,y,z]=get_coordinates_on_detector(sms.sm(k), im_par);
            if any(sms.sm(k).n_phot_ch1)
                [sms.sm(k),im_par]=update_image(x,y,z,sms.sm(k), im_par, sm_par, 1);
            end
            if im_par.two_channel==1 && any(sms.sm(k).n_phot_ch2)
                [sms.sm(k),im_par]=update_image(x,y,z,sms.sm(k), im_par, sm_par, 2);
            end
        else
            [x,y,z]=get_sub_coordinates_on_detector(sms.sm(k), im_par);
            if any(sms.sm(k).n_phot_ch1)
                [sms.sm(k),im_par]=update_sub_image(x,y,z,sms.sm(k), im_par, sm_par, 1);
            end
            if im_par.two_channel==1 && any(sms.sm(k).n_phot_ch2)
                [sms.sm(k),im_par]=update_sub_image(x,y,z,sms.sm(k), im_par, sm_par, 2);
            end
        end
    end
end
