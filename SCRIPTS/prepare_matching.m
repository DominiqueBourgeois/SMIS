function sms=prepare_matching(sms, sm_par, fluorophore_pairs)

n_pairs=size(fluorophore_pairs,1);

for i=1:n_pairs
    dye_id_ref=fluorophore_pairs(i,1); % The reference
    dye_id_target=fluorophore_pairs(i,2); % The td dyes to be placed
    n_mol_ref=sm_par(dye_id_ref).n_mol;
    n_mol_target=sm_par(dye_id_target).n_mol;
    if n_mol_target>=n_mol_ref
        newVals=num2cell(horzcat([sms(dye_id_ref).sm.id],zeros(1,n_mol_target-n_mol_ref))); 
       [sms(dye_id_target).sm.matched] = newVals{:};    
       newVals=[sms(dye_id_target).sm.id]; newVals=newVals(1:n_mol_ref);
       newVals=num2cell(newVals);
       [sms(dye_id_ref).sm.matched] = newVals{:};

    else
        newVals=[sms(dye_id_ref).sm.id]; newVals=newVals(1:n_mol_target);
        newVals=num2cell(newVals); 
       [sms(dye_id_target).sm.matched] = newVals{:}; 
       
       newVals=num2cell(horzcat([sms(dye_id_target).sm.id],zeros(1,n_mol_ref-n_mol_target)));
       [sms(dye_id_ref).sm.matched] = newVals{:};
    end

end