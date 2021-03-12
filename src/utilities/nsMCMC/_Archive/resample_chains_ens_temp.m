function mcmcInfo = resample_chains_ens_temp(mcmcInfo)

% extract parameters
A_log = log(repmat(mcmcInfo.A_curr,1,1,mcmcInfo.n_chains));
pi0_log = log(permute(repmat(mcmcInfo.pi0_curr',mcmcInfo.n_chains,1),[2 1]));
nStates = mcmcInfo.nStates;
n_traces = mcmcInfo.n_traces;
n_chains = mcmcInfo.n_chains;
seq_length = mcmcInfo.seq_length;

% generate reference vector
chain_id_ref = 0:n_chains-1;
row_ref = (1:nStates)';

% determine order of sampling    
n_reps = 2;%max([2 round(200/(mcmcInfo.step-1)^2)]);    
sample_indices = randsample(repelem(1:seq_length,n_reps),n_reps*seq_length,false);

% iterate through indices to sample
for ind = sample_indices
    mcmcInfo.ind = ind;
    %%% calculate forward-backward state probabilities
    prev_probs_log = pi0_log;                  
    if ind > 1
        % calculate linear indices 
        prev_state_array = mcmcInfo.sample_chains(ind-1,:,:);
        row_col_array = (prev_state_array-1)*nStates+row_ref;
        lin_index_array = row_col_array + chain_id_ref*nStates^2;
        % extract probabilities
        prev_probs_log = A_log(lin_index_array);            
    end
    post_probs_log = pi0_log;        
    if ind < seq_length
        % calculate linear indices 
        post_state_array = mcmcInfo.sample_chains(ind+1,:,:);
        row_col_array = (row_ref-1)*nStates+post_state_array;
        lin_index_array = row_col_array + chain_id_ref*nStates^2;
        % extract probabilities
        post_probs_log = A_log(lin_index_array);             
    end
    % combine
    logL_tr = prev_probs_log + post_probs_log;

    %%% calculate fluorescence probability component

    % calculate fluo error term          
    logL_fluo = calculate_fluo_logL_ens_temp(mcmcInfo);
    
    %%% put everything together
    total_log_likelihoods = logL_fluo + logL_tr  - logsumexp(logL_fluo + logL_tr,1);        
    total_likelihoods = exp(total_log_likelihoods);    

    %%% draw new samples
    option_array = cumsum(total_likelihoods);
    rand_array = repmat(rand(1,n_chains,n_traces),nStates,1);
    mcmcInfo.sample_chains(ind,:,:) = sum(rand_array > option_array) + 1;

end   