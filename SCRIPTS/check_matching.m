function     sms=check_matching(sms, fluorophore_pairs)

% NAME:
%	Check matching between pairs of molecules level
%   Molecules that have been matched but have lost their matching partner
%   will be assigned a matched field = -1
% MODIFICATION HISTORY:
%	D.Bourgeois, April 2011.
%	D.Bourgeois, May 2015: added z coordinate for 3D.
%	D.Bourgeois, June 2019

%now check matching: some molecules will not be matched any more 
n_pairs=size(fluorophore_pairs,1);

for i=1:n_pairs
    dye_id_ref=fluorophore_pairs(i,1); % The reference
    dye_id_target=fluorophore_pairs(i,2); % The td dyes to be placed
    id=[sms(dye_id_ref).sm.id];
    matched=[sms(dye_id_target).sm.matched];
    %search in match array those elements that are not in id   
    matched(~ismember(matched,id) & matched>0)=-1;
    newVals = num2cell(matched); [sms(dye_id_target).sm.matched] = newVals{:};
    disp(['Number of matched molecules for dye ',num2str(dye_id_target), ' is: ',num2str(nnz(matched))]);

    
    %repeat the other way around
    id=[sms(dye_id_target).sm.id];
    matched=[sms(dye_id_ref).sm.matched];
    matched(~ismember(matched,id) & matched>0)=-1;
    newVals = num2cell(matched); [sms(dye_id_ref).sm.matched] = newVals{:};
    disp(['Number of matched molecules for dye ',num2str(dye_id_ref), ' is: ',num2str(nnz(matched))]);

end

end