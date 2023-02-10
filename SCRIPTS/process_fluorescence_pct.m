function sm_cell = process_fluorescence_pct(sm_cell, sm_par, im_par)

%
% PURPOSE:
%   Main script to process flurescence emitted from single molecules
%
% INPUTS:
%   sm_cell: the single molecules (with coordinates on high-resolution image) in raster units
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   sm_cell: the single molecules updated for photophysical state
%
% MODIFICATION HISTORY:
%	D.Bourgeois, October 2019.
%	D.Bourgeois, July 2020.
%	D.Bourgeois, June 2022, adapted to parallel computing
%	D.Bourgeois, September 2022, optimization of parallel computing
%	D.Bourgeois, February 2023, remove use of parallel computing for this routine, found to be typically counterproductive in most cases.

%Extract the needed fields for sm_cell
sm=sm_cell(sm_par.idx.process_fluorescence.indices,:);

% Each mol will be processed individually
if im_par.fret_on==0 % Case of no FRET
    %parfor (k=1:sm_par.n_mol_eff,im_par.parforArg) % parfor not efficient
    %here
    for k=1:sm_par.n_mol_eff
        sm(:,k)=get_fluorescence_pct(sm(:,k), sm_par, im_par);
    end
else % Case of FRET
    %parfor (k=1:sm_par.n_mol_eff, im_par.parforArg) % parfor not efficient
    %here
    for k=1:sm_par.n_mol_eff
        sm(:,k)=get_fluorescence_in_fret_mode_pct(sm(:,k), sm_par, im_par);
    end
end

%Fill up sm_cell
sm_cell(sm_par.idx.process_fluorescence.indices,:)=sm;

