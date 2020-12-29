function mcmcInfo = resample_chains_iid(mcmcInfo)

% Assumes observations are iid to get initial estimates for V and sigma

% extract parameters
nStates = mcmcInfo.nStates;
n_traces = mcmcInfo.n_traces;
n_chains = mcmcInfo.n_chains;
seq_length = mcmcInfo.seq_length;


% initialize
sample_chains_curr = mcmcInfo.sample_chains;

% interate through traces
for n = 1:n_traces    
    mcmcInfo.trace_id = n;
  
    % slice sampling array
    mcmcInfo.sample_chains_slice = sample_chains_curr(:,:,n);    
    
    % determine order of sampling
    sample_indices = randsample(1:seq_length,seq_length,false);
    
    % iterate through indices to sample
    for ind = sample_indices
        mcmcInfo.ind = ind;        
        
        %%% calculate fluorescence probability component
           
        % calculate fluo error term      
        logL_fluo = calculate_fluo_logL(mcmcInfo);
        
        %%% put everything together
        total_log_likelihoods = exp(logL_fluo - logsumexp(logL_fluo,1));
        
        %%% draw new samples
        option_array = cumsum(total_log_likelihoods);
        rand_array = repmat(rand(1,n_chains),nStates,1);
        mcmcInfo.sample_chains_slice(ind,:) = sum(rand_array > option_array) + 1;
                
    end
    sample_chains_curr(:,:,n) = mcmcInfo.sample_chains_slice;
end   
mcmcInfo.sample_chains = sample_chains_curr;