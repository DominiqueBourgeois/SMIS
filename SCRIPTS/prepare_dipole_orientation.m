function sms=prepare_dipole_orientation(sms,n_dyes,sm_par)

for i=1:n_dyes
    if sm_par(i).dipole_orientation==0 && sm_par(i).anisotropy==1 % if randomly oriented molecules
        sms(i).sm=orient_molecules(sms(i).sm, sm_par(i).n_mol_eff);
    end
end
end