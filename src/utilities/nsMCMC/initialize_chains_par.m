function mcmcInfo = initialize_chains_par(mcmcInfo)

    % extract parameters
    A_curr = mcmcInfo.A_curr;
    nStates = size(A_curr,1);
    pi0_curr = mcmcInfo.pi0_curr;
    seq_length = mcmcInfo.seq_length;
    n_chains = mcmcInfo.n_chains;
    n_traces = mcmcInfo.n_traces;    

    % initialize
    mcmcInfo.sample_chains = NaN(seq_length,n_chains,n_traces);

    % draw initial state
    pi0_array = cumsum(pi0_curr');
    row_index_vec = repmat(1:nStates,1,n_chains);
    depth_index_vec = repelem(1:n_chains,nStates);
    for n = 1:n_traces        
        rand_array = repmat(rand(1,n_chains),nStates,1);
        mcmcInfo.sample_chains(1,:,n) = sum(rand_array > pi0_array) + 1;%randsample(state_options,n_chains,true,pi0_curr);    
        % simulate remaining time step
        for t = 2:seq_length
            % get indices
            state_vec = repelem(mcmcInfo.sample_chains(t-1,:,n),nStates);
            linear_indices = sub2ind(size(A_curr),row_index_vec,state_vec,depth_index_vec);
            
            % simulate transitions
            A_array = cumsum(reshape(A_curr(linear_indices),nStates,n_chains));
            rand_array = repmat(rand(1,n_chains),nStates,1);
            mcmcInfo.sample_chains(t,:,n) = sum(rand_array > A_array) + 1;
        end    
    end   
