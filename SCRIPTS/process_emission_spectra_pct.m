function sm_cell=process_emission_spectra_pct(sm_cell, sm_par, im_par)

%
% PURPOSE:
%   Get emission spectra and detected number of photons 
%
% INPUTS:
%   sms: the single molecules 
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   sms: the single molecules updated for emission spectra and detected
%   photons
%
% MODIFICATION HISTORY:
%	D.Bourgeois, September 2019.
%	D.Bourgeois, September 2022, optimization of parallel computing

%generate the emission spectrum for the SM

%Extract the needed fields for sm_cell
sm=sm_cell(sm_par.idx.process_emission_spectra.indices,:);

parfor (k=1:sm_par.n_mol_eff, im_par.parforArg)
%for k=1:sm_par.n_mol_eff
    sm(:,k)=get_emission_spectra_pct(sm(:,k), sm_par, im_par);
end

%Fill up sm_cell
sm_cell(sm_par.idx.process_emission_spectra.indices,:)=sm;


