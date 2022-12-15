function [n_fret]=get_fret_statistics(sm, n_mol, generate_data_set, frame, sm_par, ...
    sm_par_2, im_par)

% NAME:
%	get_fret_statistics
% PURPOSE:
%
% INPUTS:
%
%
% OUTPUTS:
%
% MODIFICATION HISTORY:
%	D.Bourgeois, December 2011.


mol1_bleached = cell2mat({sm.bleached});
mol1_activated = cell2mat({sm.activated});

%indices of the green molecules not bleached at start of frame...
ok_green = find(~mol1_activated & (~mol1_bleached | mol1_bleached==frame));
%indices of the green molecules not bleached at start of frame.
ok_red = find(mol1_activated & (~mol1_bleached | mol1_bleached==frame));

%green (non activated) molecules
if ~isempty(ok_green) % if some green molecules in current frame
    n1_green = length(ok_green); % number of green molecules still OK at start of frame
    disp('******* Green molecules ********* ');
    disp(['Number of green molecules (Dye1): ', num2str(n1_green)]);
    disp(['Fraction of green molecules (Dye1): ', num2str(n1_green/n_mol)]);
    mol1_bl = mol1_bleached(ok_green);
    n1_bleached = length(find(mol1_bl==frame)); % molecules that bleached in current frame
    mol1_n_phot = mean(cell2mat({sm(ok_green).n_phot})); % mean # of emitted photons
    mol1_n_phot_ratio = mean(cell2mat({sm(ok_green).n_phot_ratio})); % mean fraction of photons emitted in ch1
    disp(['Number of molecules that bleached in this frame(Dye1): ', num2str(n1_bleached)]);
    disp(['Fraction of molecules that bleached in this frame(Dye1): ', num2str(n1_bleached/max(n_mol,1))]);
    disp(['Mean # of emitted photons up to this frame (reaching the detector) (Dye1): ', ...
        num2str(mol1_n_phot)]);
    disp(['Fraction # of emitted photons up to this frame reaching channel 1 (Dye1): ', ...
        num2str(mol1_n_phot_ratio)]);
    if im_par.fret_on
        mol1_fret=cell2mat({sm(ok_green).fret_eff});
        mean_fret_eff=mean(mol1_fret);
        n_fret=length(find(mol1_fret > 0));
        if n_fret > 0
            mean_fret_eff2=mean(mol1_fret(mol1_fret > 0));
        else
            mean_fret_eff2=0;
        end
        disp(['Mean FRET efficiency for all molecules (Dye1 -> Dye2): ', num2str(mean_fret_eff)]);
        disp(['Mean FRET efficiency for fretting molecules(Dye1 -> Dye2): ', num2str(mean_fret_eff2)]);
        disp(['Number of green molecules FRETTING (Dye1 -> Dye2): ', num2str(n_fret)]);
        disp(['Fraction of green molecules FRETTING (Dye1 -> Dye2): ', num2str(n_fret/max(n_mol,1))]);
    else
        n_fret=0; 
    end
else % no green molecule in current frame
    n_fret=0;
    disp('No green molecule !');
end

%red (activated) molecules
if sm_par.photoactivatable
    if ~isempty(ok_red) % if some red molecules in current frame
        n1_red = length(ok_red); % number of red molecules still OK at start of frame
        disp('******* Red molecules ********* ');
        disp(['Number of red molecules (Dye1): ', num2str(n1_red)]);
        disp(['Fraction of red molecules (Dye1): ', num2str(n1_red/n_mol)]);
        mol1_bl = mol1_bleached(ok_red);
        n1_bleached = length(find(mol1_bl==frame)); % molecules that bleached in current frame
        mol1_n_phot = mean(cell2mat({sm(ok_red).n_phot})); % mean # of emitted photons
        mol1_n_phot_ratio = mean(cell2mat({sm(ok_red).n_phot_ratio})); % mean fraction of photons emitted in ch1
        disp(['Number of molecules that bleached in this frame(Dye1): ', num2str(n1_bleached)]);
        disp(['Fraction of molecules that bleached in this frame(Dye1): ', num2str(n1_bleached/max(n_mol,1))]);
        disp(['Mean # of emitted photons up to this frame (reaching the detector) (Dye1): ', ...
            num2str(mol1_n_phot)]);
        disp(['Fraction # of emitted photons up to this frame reaching channel 1 (Dye1): ', ...
            num2str(mol1_n_phot_ratio)]);
        % check # of switched red dyes       
        if sm_par_2.photoactivatable
            mol1_switched = cell2mat({sm.switched});
            n1_switched = length(find(mol1_switched==frame)); % molecules that switched in current frame
            disp(['Number of molecules that switched in this frame(Dye1): ', num2str(n1_switched)]);
            disp(['Fraction of molecules that switched in this frame(Dye1): ', num2str(n1_switched/max(n_mol,1))]);
        end
    else % no red molecule in current frame
        disp('No red molecule  !');
    end
end

if generate_data_set==1
    tot_bl = length(find(mol1_bleached >= 1));
    if sm_par.photoactivatable
        tot_act = length(find(mol1_activated >= 1));
        disp(['Total fraction of molecules activated (Dye1): ', num2str(tot_act/max(n_mol,1))]);
    end
    disp(['Total fraction of molecules bleached (Dye1): ', num2str(tot_bl/max(n_mol,1))]);
end
end




