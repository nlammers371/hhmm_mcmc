function mcmcInfo = initialize_chains(mcmcInfo)

    % extract parameters
    A_curr = mcmcInfo.A_curr;
    nStates = size(A_curr,1);
    pi0 = mcmcInfo.pi0_curr;
    seq_length = mcmcInfo.seq_length;
    n_chains = mcmcInfo.n_chains;
    n_traces = mcmcInfo.n_traces;
    state_options = 1:nStates;

    % initialize
    mcmcInfo.sample_chains = NaN(seq_length,n_chains,n_traces);

    % draw initial state
    for n = 1:n_traces
        mcmcInfo.sample_chains(1,:,n) = randsample(state_options,n_chains,true,pi0);    
        % simulate remaining time step
        for t = 2:seq_length
            % simulate transitions
            A_array = cumsum(A_curr(:,mcmcInfo.sample_chains(t-1,:,n)));
            rand_array = repmat(rand(1,n_chains),nStates,1);
            mcmcInfo.sample_chains(t,:,n) = sum(rand_array > A_array) + 1;
        end    
    end   
