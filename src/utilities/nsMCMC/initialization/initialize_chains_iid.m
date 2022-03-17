function mcmcInfo = initialize_chains_iid(mcmcInfo)
    
    % extract parameters    
    nStates = mcmcInfo.nStates;
    pi0 = mcmcInfo.pi0_curr;
    sql = mcmcInfo.seq_length;
    n_chains = mcmcInfo.n_chains_eff;
    n_traces = mcmcInfo.n_traces;    
    seq_length = mcmcInfo.seq_length;    
    
    % initialize
    mcmcInfo.sample_chains = NaN(sql, n_chains, n_traces);

    % draw from initial state PDF 
    for n = 1:n_chains
        for t = 1:n_traces
            mcmcInfo.sample_chains(:,n,t) = randsample(1:nStates,seq_length,true,pi0(n,:));
        end
    end