function sms = process_fluorescence(sms, sm_par, im_par)

%
% PURPOSE:
%   Main script to process flurescence emitted from single molecules
%
% INPUTS:
%   sms: the single molecules (with coordinates on high-resolution image) in raster units
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   sms: the single molecules updated for photophysical state
%
% MODIFICATION HISTORY:
%	D.Bourgeois, October 2019.
%	D.Bourgeois, July 2020.
%	D.Bourgeois, June 2022, adapted to parallel computing

sm=sms.sm;

% Each mol will be processed individually
if im_par.fret_on==0 % Case of no FRET
    parfor (k=1:sm_par.n_mol_eff,im_par.parforArg)
        sm(k)=get_fluorescence(sm(k), sm_par, im_par);
    end
else % Case of FRET
    parfor (k=1:sm_par.n_mol_eff, im_par.parforArg)
        sm(k)=get_fluorescence_in_fret_mode(sm(k), sm_par, im_par);
    end
end

sms.sm=sm;


