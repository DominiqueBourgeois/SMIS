function smis_par=check_smis_versions(smis_par)

%
% PURPOSE:
%   Ensure compatibility of the SMIS parameters with older SMIS versions
%
% INPUTS:
%   smis_par: the loaded SMIS parameters possibly from an older version of SMIS
%
% OUTPUTS:
%   smis_par: The updated SMIS parameters
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2024


% For SMIS 2.3

if ~isfield(smis_par,'in_images')
    smis_par.in_images=repmat({[]},smis_par.n_fluorophores,1);
end

if ~isfield(smis_par.BG,'textured_bg_pattern')
    smis_par.BG.textured_bg_pattern=[];
end

if ~isfield(smis_par.lasers,'mask_pattern')
    [smis_par.lasers.('mask_pattern')]=deal([]);
end

for k=1:smis_par.n_fluorophores
    if ~isfield(smis_par.Fluorophores(k).Motion,'Hop_Diffusion')
        smis_par.Fluorophores(k).Motion.Hop_Diffusion=0;
        smis_par.Fluorophores(k).Motion.Hop_Probability=zeros(1,numel(smis_par.Fluorophores(k).Motion.D));
    end
end

if ~isfield(smis_par.lasers,'beam_profile_dir')
    [smis_par.lasers.('beam_profile_dir')]=deal('');
    [smis_par.lasers.('beam_profile_file')]=deal('');
end


% Update drift smis vsn2.3
if ~isfield(smis_par.drift,'x1')
    smis_par.drift.x1=smis_par.drift.x_drift(1);
    smis_par.drift.x2=smis_par.drift.x_drift(2);
    smis_par.drift.x3=smis_par.drift.x_drift(3);
    smis_par.drift.x_n=smis_par.drift.x_drift(4);
    smis_par.drift.x_eval='';

    smis_par.drift.y1=smis_par.drift.y_drift(1);
    smis_par.drift.y2=smis_par.drift.y_drift(2);
    smis_par.drift.y3=smis_par.drift.y_drift(3);
    smis_par.drift.y_n=smis_par.drift.y_drift(4);
    smis_par.drift.y_eval='';

    smis_par.drift.z1=smis_par.drift.z_drift(1);
    smis_par.drift.z2=smis_par.drift.z_drift(2);
    smis_par.drift.z3=smis_par.drift.z_drift(3);
    smis_par.drift.z_n=smis_par.drift.z_drift(4);
    smis_par.drift.z_eval='';

    smis_par.drift.rot_x0=smis_par.drift.rot_drift(1);
    smis_par.drift.rot_y0=smis_par.drift.rot_drift(2);
    smis_par.drift.rot_theta=smis_par.drift.rot_drift(3);
    smis_par.drift.rot_n=smis_par.drift.rot_drift(4);  
    smis_par.drift.rot_eval='';

    smis_par.drift.dx=[];
    smis_par.drift.dy=[];
    smis_par.drift.dz=[];
    smis_par.drift.dtheta=[];

    smis_par.drift=rmfield(smis_par.drift,'x_drift');
    smis_par.drift=rmfield(smis_par.drift,'y_drift');
    smis_par.drift=rmfield(smis_par.drift,'z_drift');
    smis_par.drift=rmfield(smis_par.drift,'rot_drift');
end

% Update z_drift_range for 3D PSF calculation smis vsn2.3
if ~isfield(smis_par.drift,'z_drift_range')
    if smis_par.drift.state==1
        cz=cumsum(smis_par.drift.dz); % Cumulative drift in z
        smis_par.drift.z_drift_range=[min(cz),max(cz)];
    else
        smis_par.drift.z_drift_range=[0,0];
    end
end

% Update sample height corresponding to coverslip position smis vsn2.3
% app.ObjectivePSFpar.sample_z_coverslip=value;
if ~isfield(smis_par.obj_and_psf,'sample_z_coverslip')
    smis_par.obj_and_psf.sample_z_coverslip=0;
end

