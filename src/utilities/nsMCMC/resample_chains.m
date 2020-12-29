function mcmcInfo = resample_chains(mcmcInfo)

% extract parameters
A_curr = mcmcInfo.A_curr;
A_log = log(mcmcInfo.A_curr);
nStates = size(A_curr,1);
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
        %%% calculate forward-backward state probabilities
        prev_probs_log = 0;                  
        if ind > 1
            prev_probs_log = A_log(:,mcmcInfo.sample_chains_slice(ind-1,:));            
        end
        post_probs_log = 0;        
        if ind < seq_length
            post_probs_log = A_log(mcmcInfo.sample_chains_slice(ind+1,:),:);            
        end
        % combine
        logL_tr = prev_probs_log + post_probs_log';
        
        %%% calculate fluorescence probability component
           
        % calculate fluo error term      
        logL_fluo = calculate_fluo_logL(mcmcInfo);
        
        %%% put everything together
        total_log_likelihoods = logL_fluo + logL_tr;
        total_log_likelihoods = exp(total_log_likelihoods - logsumexp(total_log_likelihoods,1));
        
        %%% draw new samples
        option_array = cumsum(total_log_likelihoods);
        rand_array = repmat(rand(1,n_chains),nStates,1);
        mcmcInfo.sample_chains_slice(ind,:) = sum(rand_array > option_array) + 1;
                
    end
    sample_chains_curr(:,:,n) = mcmcInfo.sample_chains_slice;
end   
mcmcInfo.sample_chains = sample_chains_curr;