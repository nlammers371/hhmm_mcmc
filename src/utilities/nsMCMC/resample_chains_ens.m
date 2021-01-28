function mcmcInfo = resample_chains_ens(mcmcInfo)

% extract parameters
A_log = log(mcmcInfo.A_curr);
pi0 = mcmcInfo.pi0_curr;
pi0_log = log(pi0);
nStates = mcmcInfo.nStates;
n_traces = mcmcInfo.n_traces;
n_chains = mcmcInfo.n_chains;
nSteps = mcmcInfo.nSteps;
seq_length = mcmcInfo.seq_length;

% generate reference vector
mcmcInfo.chain_id_ref = 0:n_chains-1;
mcmcInfo.trace_id_ref = reshape(0:n_traces-1,1,1,[]);
mcmcInfo.row_ref = (1:nStates)';
mcmcInfo.step_ref = (-nSteps+1:nSteps-1)';

% generate random sampling orders
n_reps = 1;
sample_indices = NaN(n_reps*seq_length,size(mcmcInfo.sample_chains,2),size(mcmcInfo.sample_chains,3));
for n = 1:n_chains
    for t = 1:n_traces
        sample_indices(:,n,t) = randsample(repelem(1:seq_length,n_reps),n_reps*seq_length,false);
    end
end

% generate temporary chain array
sample_chains_temp = NaN(size(mcmcInfo.sample_chains,1)+2,size(mcmcInfo.sample_chains,2),size(mcmcInfo.sample_chains,3));
sample_chains_temp(2:end-1,:,:) = mcmcInfo.sample_chains;

% initialize starting and end indices randomly using pi0
sample_chains_temp(1,:,:) = reshape(randsample(1:nStates,n_chains*n_traces,true,pi0),1,n_chains,n_traces);
sample_chains_temp(end,:,:) = reshape(randsample(1:nStates,n_chains*n_traces,true,pi0),1,n_chains,n_traces);
seq_len_temp = seq_length+2;
seq_len_dummy = seq_length+2*(nSteps-1);
% create dummy chain array for fluorescence sampling
% true_states = repmat(permute(vertcat(mcmcInfo.masterSimStruct.naive_states)',[1 3 2]),1, n_chains,1);

% true_states =repmat(reshape(vertcat(mcmcInfo.masterSimStruct.naive_states)',seq_length,1,n_traces),1,n_chains,1);
% warning('using actua; states for testing')
mcmcInfo.sample_chains_dummy = cat(1,ones(nSteps-1,n_chains,n_traces),mcmcInfo.sample_chains,ones(nSteps-1,n_chains,n_traces));
mcmcInfo.observed_fluo_dummy = cat(1,zeros(nSteps-1,n_traces),mcmcInfo.observed_fluo,zeros(nSteps-1,n_traces));

% iterate through indices to sample
for i = 1:seq_length
    mcmcInfo.indArray = sample_indices(i,:,:);    
    
    prev_time_index_array = sample_indices(i,:,:);    
    prev_lin_index_array = prev_time_index_array + mcmcInfo.chain_id_ref*seq_len_temp + mcmcInfo.trace_id_ref*seq_len_temp*n_chains;
    post_lin_index_array = prev_lin_index_array + 2;
    
    %%% previous state %%%
    % calculate linear indices 
    prev_state_array = sample_chains_temp(prev_lin_index_array);
    row_col_array_from = (prev_state_array-1)*nStates+mcmcInfo.row_ref;    
    
    % extract probabilities
    prev_probs_log = A_log(row_col_array_from);
          
    %%% following state %%%
    % calculate linear indices 
    post_state_array = sample_chains_temp(post_lin_index_array);  
    row_col_array_to = (mcmcInfo.row_ref-1)*nStates+post_state_array;
        
    % extract probabilities
    post_probs_log = A_log(row_col_array_to);
       
    % combine
    logL_tr = prev_probs_log + post_probs_log;

    %%% calculate fluorescence probability component

    % calculate fluo error term          
    logL_fluo = calculate_fluo_logL_ens(mcmcInfo);
    
    %%% put everything together
    total_log_likelihoods = logL_tr + logL_fluo;
    total_log_likelihoods = total_log_likelihoods - logsumexp(total_log_likelihoods,1);
    total_likelihoods = exp(total_log_likelihoods);    

    %%% draw new samples
    option_array = cumsum(total_likelihoods);
    rand_array = repmat(rand(1,n_chains,n_traces),nStates,1);
%     
%     mcmcInfo.sample_chains(curr_lin_index_array) = sum(rand_array > option_array) + 1;
    sample_chains_temp(prev_lin_index_array+1) = sum(rand_array > option_array) + 1;
    curr_lin_index_array = (mcmcInfo.indArray+nSteps-1) + mcmcInfo.chain_id_ref*seq_len_dummy + mcmcInfo.trace_id_ref*seq_len_dummy*n_chains;
    mcmcInfo.sample_chains_dummy(curr_lin_index_array) = sample_chains_temp(prev_lin_index_array+1);

end  

mcmcInfo.sample_chains = sample_chains_temp(2:end-1,:,:);