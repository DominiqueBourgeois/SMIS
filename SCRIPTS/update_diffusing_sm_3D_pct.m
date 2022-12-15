function  sm_cell = update_diffusing_sm_3D_pct(sm_cell, im_par, sm_par, S_DS)

% PURPOSE:
%	Update 3D position and diffusion state of diffusing SM that are in
%	starting diffusion state S_DS
%
% INPUTS:
%   sm_cell: the single molecules 
%	sm_par: the sm parameters
%	im_par: the imaging parameters
%   S_DS: starting diffusion state
%
% OUTPUTS:
%   sm_cell: the single molecules updated for: diffusion states, current position
%
% MODIFICATION HISTORY:
%	D.Bourgeois, December 2019: version > simulate_palm_vsn15.3
%	D.Bourgeois, November 2020: version > simulate_palm_vsn16.3 (introduce
%	Velocity)
%	D.Bourgeois, January 2022: modify script to account for fast exchange 
%       rates between diffusion states in the diffuse PSF mode. Adapted to SMIS 1.3
%	D.Bourgeois, April 2022: option to record the whole diffusion state
%       history, but not used at this stage.
%	D.Bourgeois, September 2022, optimized for parallel computing

%Extract the needed fields for sm_cell
sm=sm_cell(sm_par.w_idx,:);

% Do not record the ds history. May be useful for a future version
sm_par.record_ds_history=0; % 

% Check if exchange rates are compatible with frame time and add time
if im_par.use_diffuse_psf==1
    if ~isempty(sm_par.D_ex_rates)
        max_ex_rate=max(sm_par.D_ex_rates(:,3));
        %Get the maximum subframe time [s] with oversampling of ex_rates_min_oversampling
        sm_par.max_dt=1/(max_ex_rate*im_par.ex_rates_min_oversampling);
    else
        sm_par.max_dt=inf;
    end
else
    sm_par.max_dt=[]; % Needed for parallel computing
end


n_mol=size(sm_cell,2);
parfor (k=1:n_mol, im_par.parforArg) % Go for all activated molecules
%for k=1:n_mol
    sm(:,k)=move_one_diffusing_sm_3D_pct(sm(:,k), sm_par, im_par, S_DS);
end

%Fill up sm_cell
sm_cell(sm_par.w_idx,:)=sm;



