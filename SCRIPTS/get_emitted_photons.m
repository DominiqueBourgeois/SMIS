function n = get_emitted_photons(n_tot, theta, im_par, sm_par)
% NAME:
%   get_emitted_photons
%
% PURPOSE:
%   Calculates #of photons reaching objective based on Fourkas, Opt. Lett. 2001
%
% INPUTS:
%   n_tot: total # of emitted photons
%   theta: polar angle of molecule [degrees]
%   im_par: imaging parameters
%
% OUTPUTS:
%	a: efficiency of collection
%
% MODIFICATION HISTORY:
%	D.Bourgeois, May 2012. September 2019
%

%Only a fraction will reach detector
%Also take into account here detector quantum efficiency
if im_par.obj.eff_from_opening_angle==1 && sm_par.anisotropy==1
    %  Furkas formula
    n=fix(im_par.det.QE*2*im_par.obj.mic_transmission*n_tot*(im_par.obj.eff_a+im_par.obj.eff_b*(cos(theta*pi/180))^2));
else
    n=fix(im_par.det.QE*n_tot*im_par.obj.mean_det_eff);
end
if im_par.debug
    disp(['Detected photons: ', num2str(n)]);
end
end