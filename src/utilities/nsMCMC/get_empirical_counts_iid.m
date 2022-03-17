function mcmcInfo = get_empirical_counts_iid(mcmcInfo)
    
    % extract parameters
    pi0_curr = mcmcInfo.pi0_curr;
    nStates = size(pi0_curr,2);
    n_chains = mcmcInfo.n_chains_eff;


    % initialize count arrays     
    mcmcInfo.state_count_array = zeros(n_chains,nStates);       
    
    % get state counts
    for k = 1:nStates
        mcmcInfo.state_count_array(:,k) = sum(sum(mcmcInfo.sample_chains==k,3),1);
    end        
    