function C=prepare_XYZ(n_dyes, sms, im_par)

% PURPOSE:
%	Prepare XYZ coordinates from sms
%
% INPUTS:
%	sms: the sms
%   im_par: image parameters
%
% OUTPUTS:
%   C = the X, Y and Z coordinates; xd, yd and zd are also defined: they
%   will contain the moving coordinates along data set in case of drift or
%   diffusion
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15

C(1:n_dyes)=struct('x', [], 'y', [], 'z', []);
for i=1:n_dyes
    %Get the xy coordinates on the binned image
    x=[sms(i).sm.x]; y=[sms(i).sm.y]; z=[sms(i).sm.z];
    % for example if binning=4 and x goes from 0.5 to 512.5 in unbinned image, it has to
    % go from 0.5 to 128.5 in binned image
    C(i).x=(x-0.5)/im_par.binning+0.5;
    C(i).y=(y-0.5)/im_par.binning+0.5;
    
    if im_par.simul_3D==1
        C(i).z=(z-0.5)/im_par.binning+0.5;
    end   
end
end