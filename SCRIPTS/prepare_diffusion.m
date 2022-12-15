function sm_par=prepare_diffusion(n_dyes, sm_par)

for i=1:n_dyes
    sm_par(i).D_rate_matrix=get_diff_rates_matrix(numel(sm_par(i).D),sm_par(i).D_ex_rates);
end


