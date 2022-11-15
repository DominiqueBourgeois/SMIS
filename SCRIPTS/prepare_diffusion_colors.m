function  diffusion_colors=prepare_diffusion_colors(n_dyes, sm_par)

% PURPOSE:
%	Prepare diffusion image
%
% INPUTS:
%	sm_par: the sm parameters
%
% OUTPUTS:
%   diffusion_colors = random rgb diffusion colors for each molecule
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15

diffusion_colors(1:n_dyes)=struct('c',[]);
for i=1:n_dyes
    diffusion_colors(i).c=rand(sm_par(i).n_mol_eff,3); % pick random set of colors
end
end