function [sms,a_sm_all,sm_par,ok] = place_single_molecules(sms, a_all,sm_par,im_par)
% NAME:
%	PLACE_SINGLE_MOLECULES
%
% PURPOSE:
%	Get xyz coordinates of SMs following a pattern on an image
%
% CATEGORY:
%	Single molecule.
%
% CALLING SEQUENCE:
%	[sms,a_sm_all] = place_single_molecules(sms, a_all,sm_par,im_par)
%
% INPUTS:
%	sms: the single molecules
%	a_all: the images onto place the SMs
%	sm_par: general parameters on the SMs
%	im_par: general parameters on the imaging experiment
%
% OUTPUTS:
%	sms = the modified single molecules
%   a_sm_all = the image, where the SMs have been positionned
%	sm_par: updated general parameters on the SMs
%   ok = 1 if execution was okay, 0 if an error occurred
%
% COMMON BLOCKS:
%	None.
%
% SIDE EFFECTS:
%	None.
%
% RESTRICTIONS:
%	None.
%
% MODIFICATION HISTORY:
%	D.Bourgeois, April 2011.
%	D.Bourgeois, May 2015: added z coordinate for 3D.
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15
%	D.Bourgeois, March 2021: Added error handling

n_dyes=size(a_all,1);
a_sm_all=a_all;

for i=1:n_dyes
    if sm_par(i).is_acceptor==0 % Only place molecules that are not tandem targets (these will be placed next in routine place_single_td_molecules)
        if  im_par.simul_3D==0  % 2D image       
            [x_h,y_h,~,sp, a_sm_all(i,:,:),sm_par(i),ok] = place_single_molecules_simple(squeeze(a_all(i,:,:)),im_par,sm_par(i));           
            if ~ok; return; end
            newVals=num2cell(x_h); [sms(i).sm.x]=newVals{:};
            newVals=num2cell(y_h); [sms(i).sm.y]=newVals{:};
        else
            [x_h,y_h,z_h,sp, a_sm_all(i,:,:,:),sm_par(i),ok] = place_single_molecules_simple(squeeze(a_all(i,:,:,:)),im_par,sm_par(i));
            if ~ok; return; end
            newVals=num2cell(x_h); [sms(i).sm.x]=newVals{:};
            newVals=num2cell(y_h); [sms(i).sm.y]=newVals{:};
            newVals=num2cell(z_h); [sms(i).sm.z]=newVals{:};
        end 
        newVals=num2cell(sp); [sms(i).sm.c_sp]=newVals{:};
    end
end
end