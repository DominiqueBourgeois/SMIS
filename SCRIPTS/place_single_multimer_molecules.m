function     [sms,a_sm_all,ok]=place_single_multimer_molecules(sms, a_all, sm_par, im_par, n_clusters)

% NAME:
%	PLACE_SINGLE_MOLECULES
%
% PURPOSE:
%	Get xy coordinates of SMs following a pattern on an image
%
% CATEGORY:
%	Single molecule.
%
% INPUTS:
%	a_all: the image onto place the SMs
%	sm_par: sm parameters
%   n_clusters: # of clusters
%
% OUTPUTS:
%	sms = the updated sms 
%   a_sm_all = the image, where the SMs have been positionned
%   ok = 1 if execution was okay, 0 if an error occurred
%
% MODIFICATION HISTORY:
%	D.Bourgeois, April 2011.
%	D.Bourgeois, May 2015: added z coordinate for 3D.
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15
%	D.Bourgeois, March 2021: Added error handling

n_dyes=size(a_all,1);
% C(1:n_dyes)=struct('x',[],'y',[],'z',[]);
% a_sm_all=zeros(size(a_all));
a_sm_all=a_all;

for i=1:n_dyes
    if sm_par(i).is_acceptor==0 % Only place molecules that are not tandem targets (these will be placed next)
        if  im_par.simul_3D==0  % 2D image
            [X,Y,~,a_sm_all(i,:,:),ok] = place_single_multimer_molecules_simple(squeeze(a_all(i,:,:)),im_par,sm_par(i), n_clusters);
            if ~ok; return; end
            newVals=num2cell(X); [sms(i).sm.x]=newVals{:};
            newVals=num2cell(Y); [sms(i).sm.y]=newVals{:};
        else
            [X,Y,Z,a_sm_all(i,:,:,:),ok] = place_single_multimer_molecules_simple(squeeze(a_all(i,:,:,:)),im_par,sm_par(i), n_clusters);
            if ~ok; return; end
            newVals=num2cell(X); [sms(i).sm.x]=newVals{:};
            newVals=num2cell(Y); [sms(i).sm.y]=newVals{:};
            newVals=num2cell(Z); [sms(i).sm.z]=newVals{:};
        end  
    end
end
end