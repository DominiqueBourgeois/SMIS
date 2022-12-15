function true_3D_kernel=prepare_3D_kernel(a_h_all)

% PURPOSE:
%	Prepare a 3D kernel with all true positions of all sms
%
% INPUTS:
%	a_h_all: the images (3D kernels) onto the SMs are placed
%
% OUTPUTS:
%   true_3D_kernel = the 3D kernel image, where all SMs have been positionned
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15.3 Corrected bug:
%	true_3D_kernel cannot be 4D array

n_dyes=size(a_h_all, 1);

true_3D_kernel=squeeze(a_h_all(1,:,:,:)>0); % set to 1 for the first dye
if n_dyes>1
    for i=2:n_dyes
        true_3D_kernel=true_3D_kernel+squeeze(a_h_all(i,:,:,:)>0); % add i for the ith dye
    end
end
end