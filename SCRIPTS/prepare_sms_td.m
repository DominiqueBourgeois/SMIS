function sms_td=prepare_sms_td(n_dyes, sms, fluorophore_pairs)

%Define tandem single molecules
sms_td(1:n_dyes)=struct('sm', []);
for i=1:n_dyes  
    sms_td(i).sm=sms(fluorophore_pairs(i)).sm;
end