function mcmcInfo = get_empirical_counts_PF(mcmcInfo)
    
    % transition counts  
    pss_array = mcmcInfo.pss_array;
    mcmcInfo.transition_count_array = sum(sum(pss_array,3),4);    
    mcmcInfo.pi0_counts = NaN(mcmcInfo.nStates,1);
    
    % initial state array
    first_states = mcmcInfo.sample_chains(1,:,:);
    first_probs = exp(mcmcInfo.sample_logL_bkd(1,:,:));
    for k = 1:mcmcInfo.nStates
        mcmcInfo.pi0_counts(k) = sum(first_probs(first_states==k)); 
    end
    