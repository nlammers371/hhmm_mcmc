function mcmcInfo = initialize_chains_ens(mcmcInfo)

    % extract parameters
    A_curr = mcmcInfo.A_curr;
    nStates = mcmcInfo.nStates;
    pi0 = mcmcInfo.pi0_curr;
    sql = mcmcInfo.seq_length;
    n_chains = mcmcInfo.n_chains;
    n_traces = mcmcInfo.n_traces;    

    % generate reference vector
    row_ref = (1:nStates)';
    % initialize
    mcmcInfo.sample_chains = NaN(sql, n_chains, n_traces);

    % draw from initial state PDF 
    mcmcInfo.sample_chains(1,:,:) = reshape(randsample(1:nStates,n_chains*n_traces,true,pi0),1,n_chains,n_traces);
    
    % iterate step by step
    for t = 2:sql
        % generate sampling indices
        prev_state_array = mcmcInfo.sample_chains(t-1,:,:);
        row_col_array = (prev_state_array-1)*nStates+row_ref;
        lin_index_array = row_col_array;% chain_id_ref*nStates^2;
        % simulate transitions
        A_array = cumsum(A_curr(lin_index_array));
        rand_array = repmat(rand(1,n_chains,n_traces),nStates,1);
        mcmcInfo.sample_chains(t,:,:) = sum(rand_array > A_array) + 1;
    end  
    