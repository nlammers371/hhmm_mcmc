function mcmcInfo = resample_chains_v4_mh(mcmcInfo)

% This script resamples the microscopic promoter state for each extant
% chain in an asynchronous manner. 

% extract parameters
A_log = log(mcmcInfo.A_curr);
pi0 = mcmcInfo.pi0_curr;
nStates = mcmcInfo.nStates;
n_traces = mcmcInfo.n_traces;
n_chains = mcmcInfo.n_chains_eff;
nStepsMax = mcmcInfo.nStepsMax;
seq_length = mcmcInfo.seq_length;

% generate reference vector
mcmcInfo.chain_id_ref = 0:n_chains-1;
mcmcInfo.trace_id_ref = reshape(0:n_traces-1,1,1,[]);
mcmcInfo.row_ref = (1:nStates)';
mcmcInfo.step_ref = (-nStepsMax+1:nStepsMax-1)';

% generate random sampling orders
n_reps = mcmcInfo.n_reps;
sample_indices = randsample(repelem(1:seq_length,n_reps),n_reps*seq_length,false);

% generate temporary chain array that includes post and prior states to
% make resampling easier below
sample_chains_temp = NaN(size(mcmcInfo.sample_chains,1)+2,size(mcmcInfo.sample_chains,2),size(mcmcInfo.sample_chains,3));
sample_chains_temp(2:end-1,:,:) = mcmcInfo.sample_chains;

% initialize post and prior states randomly using pi0 PDF
for n = 1:n_chains
    sample_chains_temp(1,n,:) = reshape(randsample(1:nStates,n_traces,true,pi0(n,:)),1,1,n_traces);
    sample_chains_temp(end,n,:) = reshape(randsample(1:nStates,n_traces,true,pi0(n,:)),1,1,n_traces);
end    

% a initialize dummy chain variable
mcmcInfo.sample_chains_temp = sample_chains_temp;

% preallocate helper array for indexing
index_helper2 = mcmcInfo.chain_id_ref*nStates^2 + (mcmcInfo.row_ref-1)*nStates;

% iterate through indices to sample
for i = 1:seq_length*n_reps  
    
    mcmcInfo.samp_index = sample_indices(i);        
    prev_time_index = sample_indices(i);    
    post_time_index = prev_time_index + 2;
    
    %%% previous state %%%
    % calculate linear indices 
    prev_state_array = sample_chains_temp(prev_time_index,:,:);
    row_col_array_from = mcmcInfo.chain_id_ref*nStates^2 + (prev_state_array-1)*nStates+mcmcInfo.row_ref;    
    
    % extract probabilities
    prev_probs_log = A_log(row_col_array_from);
          
    %%% following state %%%
    % calculate linear indices 
    post_state_array = sample_chains_temp(post_time_index,:,:);  
    row_col_array_to = index_helper2 + post_state_array;
        
    % extract probabilities
    post_probs_log = A_log(row_col_array_to);
       
    % combine
    logL_tr = prev_probs_log + post_probs_log;
    
    % zero out current state
%     curr_inds = sample_chains_temp(prev_time_index+1,:,:) + mcmcInfo.chain_id_ref*nStates + mcmcInfo.trace_id_ref*nStates*n_chains;
%     logL_tr(curr_inds) = -Inf;
    
    % use transition probabilities as jump proposal distribution
    logL_tr_norm = exp(logL_tr - logsumexp(logL_tr,1));
    
    % randomly select new state proposals
    option_array = cumsum(logL_tr_norm);
    rand_array = repmat(rand(1,n_chains,n_traces),nStates,1);
    mcmcInfo.state_proposals = sum(rand_array > option_array,1) + 1;    
    
    %%% calculate fluorescence probability component
    % calculate fluo error term           
    logL_fluo = calculate_fluo_logL_v4_mh(mcmcInfo);              
    
    %%% put everything together
    logL_diffs = logL_fluo(2,:,:) - logL_fluo(1,:,:);    
    
    %%% draw new samples        
    move_flags = exp(logL_diffs) > rand(1,n_chains,n_traces);
    
    % assign new state selections
    curr_slice = sample_chains_temp(prev_time_index+1,:,:);
    prop_slice = mcmcInfo.state_proposals;
    curr_slice(move_flags) = prop_slice(move_flags);
    sample_chains_temp(prev_time_index+1,:,:) = curr_slice;
    
    % pass along to dummy array    
    mcmcInfo.sample_chains_temp = sample_chains_temp;
   
end 

mcmcInfo.sample_chains = sample_chains_temp(2:end-1,:,:);