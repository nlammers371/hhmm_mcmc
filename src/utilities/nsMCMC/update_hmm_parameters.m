function mcmcInfo = update_hmm_parameters(mcmcInfo)    

    % update A
    A_counts = sum(mcmcInfo.transition_count_array,3)/mcmcInfo.n_chains;    
    mcmcInfo.A_curr = sample_A_dirichlet(mcmcInfo.A_alpha, A_counts);
    mcmcInfo.A_inf_array(:,:,mcmcInfo.step) = mcmcInfo.A_curr;