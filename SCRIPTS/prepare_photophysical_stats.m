function photophysical_stats=prepare_photophysical_stats(n_dyes, sm_par, im_par)

% PURPOSE:
%	Prepare photophysical statistics from sms
%
% INPUTS:
%	n_dyes: the # of dyes
%   sm_par: sm parameters
%   im_par: image parameters
%
% OUTPUTS:
%   photophysical_stats = the photophysical statistics
%
% MODIFICATION HISTORY:
%	D.Bourgeois, June 2019: version > simulate_palm_vsn15

photophysical_stats(1:n_dyes)=struct('n_mol',[]);

for i=1:n_dyes
    n_states=sm_par(i).n_states;
    fluorescent_states=sm_par(i).fluorescent_states;
    bleached_states=sm_par(i).bleached_states;
    dark_states=find(~ismember(1:n_states,horzcat([fluorescent_states, bleached_states]))==1);  
    clear('n_mol');
    n_mol=struct(...
        'fluorescent_states', zeros(numel(fluorescent_states),im_par.n_images), ... %  # of fluorescent dye molecules in frame
        'bleached_states', zeros(numel(bleached_states),im_par.n_images), ... %  # of bleached dye molecules in frame
        'dark_states', zeros(numel(dark_states),im_par.n_images), ... %  # of dark dye molecules in frame
        'cum_fluorescent_states', zeros(numel(fluorescent_states),im_par.n_images), ... %  # of cumulated fluorescent dye molecules in frame
        'cum_bleached_states', zeros(numel(bleached_states),im_par.n_images), ... %  # of cumulated bleached dye molecules in frame
        'cum_dark_states', zeros(numel(dark_states),im_par.n_images) ... %  # of cumulated dark dye molecules in frame
        );
    photophysical_stats(i).n_mol=n_mol;
end
end