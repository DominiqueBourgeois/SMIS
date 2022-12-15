function [sms,sm_par,a_h_all] = check_maturation_level(sms, sm_par, im_par, a_h_all)
% NAME:
%	Check maturation level
%
% MODIFICATION HISTORY:
%	D.Bourgeois, April 2011.
%	D.Bourgeois, May 2015: added z coordinate for 3D.
%	D.Bourgeois, June 2019

n_dyes=size(sms,2);

for i=1:n_dyes
    sm_par(i).n_mol_eff=round(sm_par(i).maturation_level*sm_par(i).n_mol);
    if sm_par(i).maturation_level<1
        sm_select=randperm(sm_par(i).n_mol,sm_par(i).n_mol_eff);
        sms(i).sm=sms(i).sm(sm_select);
        if im_par.simul_3D==0 % 2D mode
            a_h_all(i,:,:)=0;  
            x=round([sms(i).sm.x]);
            y=round([sms(i).sm.y]);
            a_h_all(sub2ind(size(a_h_all),i*ones(1,numel(x)),x,y))=100;
        else % 3D mode
            a_h_all(i,:,:,:)=0;
            x=round([sms(i).sm.x]);
            y=round([sms(i).sm.y]);
            z=round([sms(i).sm.z]);
            a_h_all(sub2ind(size(a_h_all),i*ones(1,numel(x)),x,y,z))=100;
        end
    end
end

end