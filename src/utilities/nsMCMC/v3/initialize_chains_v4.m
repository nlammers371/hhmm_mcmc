function mcmcInfo = initialize_chains_v4(mcmcInfo)
    
    % extract parameters
    A_curr = mcmcInfo.A_curr;
    nStates = mcmcInfo.nStates;
    pi0 = mcmcInfo.pi0_curr;
%     sql = mcmcInfo.seq_length;
    n_chains = mcmcInfo.n_chains;
    n_traces = mcmcInfo.n_traces;    
    n_temps = mcmcInfo.n_temps_per_chain;
    
    % generate reference vector
    row_ref = (1:nStates)';
    
    % initialize
    mcmcInfo.sample_chains = cell(1,n_traces);
    for n = 1:n_traces
        sql = length(mcmcInfo.observed_fluo{n});
        mcmcInfo.sample_chains{n} = NaN(sql, n_chains, n_temps);
    end

    % draw from initial state PDF 
    for n = 1:n_chains
        for m = 1:n_traces
            mcmcInfo.sample_chains{m}(1,n,:) = reshape(randsample(1:nStates,n_temps,true,pi0(n,:)),1,1,n_temps);
        end
    end
       
    % iterate step by step
    for m = 1:n_traces
        sql = length(mcmcInfo.observed_fluo{m});
        % extract trace array block
        trace_array_temp = mcmcInfo.sample_chains{m};                                          
        for t = 2:sql
            % generate sampling indices
            prev_state_array = trace_array_temp(t-1,:,:);
            row_col_array = (prev_state_array-1)*nStates+row_ref;
            lin_index_array = row_col_array + mcmcInfo.chain_id_ref*nStates^2;% chain_id_ref*nStates^2;

            % simulate transitions
            A_array = cumsum(A_curr(lin_index_array));
            rand_array = repmat(rand(1,1,n_temps),nStates,1);
            trace_array_temp(t,:,:) = sum(rand_array > A_array) + 1;
        end          
        mcmcInfo.sample_chains{m} = trace_array_temp;
    end        