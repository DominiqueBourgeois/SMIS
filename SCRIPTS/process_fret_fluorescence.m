function [donor, acceptor] = process_fret_fluorescence(donor, donor_par, acceptor, acceptor_par, im_par)

%
% PURPOSE:
%   Main script to process fluorescence emitted from FRET pair
%
% INPUTS:
%   donor: the donor molecules (with coordinates on high-resolution image) in raster units
%	donor_par: the donor parameters
%   acceptor: the acceptor molecules (with coordinates on high-resolution image) in raster units
%	acceptor_par: the acceptor parameters
%	im_par: the imaging parameters
%
% OUTPUTS:
%   donor : the donor molecules updated for photophysical state
%   acceptor: the acceptor molecules updated for photophysical state
%
% MODIFICATION HISTORY:
%	D.Bourgeois, July 2020.

n_pairs=min([acceptor_par.n_mol_eff,donor_par.n_mol_eff]);

%Process the eventual left over (unpaired) molecules
if acceptor_par.n_mol_eff>donor_par.n_mol_eff
    %identify the left over acceptors
    unpaired_indices_A=find([acceptor.sm.matched]==0);
    unpaired_indices_D=[];
    paired_indices_A=find([acceptor.sm.matched]>0);
    paired_indices_D=1:donor_par.n_mol_eff;
elseif donor_par.n_mol_eff>acceptor_par.n_mol_eff
    %identify the left over donors
    unpaired_indices_D=find([donor.sm.matched]==0);
    unpaired_indices_A=[];
    paired_indices_D=find([donor.sm.matched]>0);
    paired_indices_A=1:acceptor_par.n_mol_eff;
else
    unpaired_indices_D=[];
    unpaired_indices_A=[];
    paired_indices_D=1:donor_par.n_mol_eff;
    paired_indices_A=1:acceptor_par.n_mol_eff;
end

for k=1:n_pairs
    donor.sm(paired_indices_D(k))=get_fluorescence_in_fret_mode(donor.sm(paired_indices_D(k)), donor_par, im_par);
    acceptor.sm(paired_indices_A(k))=get_fluorescence_in_fret_mode(acceptor.sm(paired_indices_A(k)), acceptor_par, im_par);
end

%process left over donors
if ~isempty(unpaired_indices_D)
    for k=1:numel(unpaired_indices_D)
        donor.sm(unpaired_indices_D(k))=get_fluorescence(donor.sm(unpaired_indices_D(k)), donor_par, im_par);
    end
end

%process left over acceptors
if ~isempty(unpaired_indices_A)
    for k=1:numel(unpaired_indices_A)
        acceptor.sm(unpaired_indices_A(k))=get_fluorescence(acceptor.sm(unpaired_indices_A(k)), acceptor_par, im_par);
    end
end

end

