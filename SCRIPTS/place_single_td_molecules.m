function [sms,sm_par,a_all_sm,ok] = place_single_td_molecules(sms, a_all_sm, a_all, fluorophore_pairs,sm_par,im_par)
% NAME:
%	PLACE_SINGLE_TANDEM_MOLECULES
%
% PURPOSE:
%	Get xyz coordinates of tandem SMs following a pattern on an image
%
% CATEGORY:
%	Single molecule.
%
% CALLING SEQUENCE:
% [sms,sm_par,a_all_sm] = place_single_td_molecules(sms, a_all_sm, a_all,fluorophore_pairs,sm_par,im_par)
%
% INPUTS:
%	sms: the single molecules
%	a_all_sm: the images onto place the SMs
%	a_all: the image patterns
%   fluorophore_pairs: the pairing of SMs
%	sm_par: general parameters on the SMs
%	im_par: general parameters on the imaging experiment
%   ok = 1 if execution was okay, 0 if an error occurred
%
% OUTPUTS:
%	sms = the modified single molecules
%   a_sm_all = the image, where the SMs have been positionned

%
% MODIFICATION HISTORY:
%	D.Bourgeois, April 2011.
%	D.Bourgeois, May 2015: added z coordinate for 3D.
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15
%	D.Bourgeois, November 2019: version > simulate_palm_vsn15.3
%	D.Bourgeois, March 2021: Added error handling

n_pairs=size(fluorophore_pairs,1);

for i=1:n_pairs
    dye_id_ref=fluorophore_pairs(i,1); % The reference
    dye_id_target=fluorophore_pairs(i,2); % The td dyes to be placed
    n_mol_ref=sm_par(dye_id_ref).n_mol;
    n_mol_target=sm_par(dye_id_target).n_mol;
    if im_par.simul_3D==0
        a_ref=squeeze(a_all(dye_id_ref, :,:));
        a_target=squeeze(a_all(dye_id_target, :,:));
    elseif im_par.simul_3D==1
        a_ref=squeeze(a_all(dye_id_ref,:,:,:));
        a_target=squeeze(a_all(dye_id_target, :,:,:));
    end
    
    [sms(dye_id_ref).sm, sms(dye_id_target).sm, sm_par(dye_id_target), ...
        a_all_sm(dye_id_target, :,:), n_mol_ref_out, n_mol_target_out,ok]=...
        place_linked_dye_in_patterns(sms(dye_id_ref).sm, sms(dye_id_target).sm, sm_par(dye_id_target),...
        a_ref, a_target, im_par, n_mol_ref, n_mol_target, i);
    if ~ok; return; end
    sm_par(dye_id_target).n_mol=n_mol_target_out;
    sm_par(dye_id_ref).n_mol=n_mol_ref_out;
    
    if n_mol_target_out==0
        errordlg('Error: No linked molecule could be positionned: increase d1d2_dist !','SMIS Error');
        uiwait
        ok=0;
        return
    end
end


end