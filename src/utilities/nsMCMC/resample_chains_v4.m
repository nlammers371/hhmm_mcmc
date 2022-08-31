function mcmcInfo = resample_chains_v4(mcmcInfo)

% This script resamples the microscopic promoter state for each extant
% chain in an asynchronous manner. 
tic
% extract parameters
A_log = log(mcmcInfo.A_curr);
pi0 = mcmcInfo.pi0_curr;
nStates = mcmcInfo.nStates;
n_traces = mcmcInfo.n_traces;
n_chains = mcmcInfo.n_chains_eff;
% nStepsMax = mcmcInfo.nStepsMax;
seq_length = mcmcInfo.seq_length;

% generate lookup table with ms2 kernel values
nStepsMax = 0;
for n = 1:n_chains
    nStepsMax = max([nStepsMax find(mcmcInfo.coeff_MS2(:,n)>0,1,'last')]);
end    
mcmcInfo.max_w = nStepsMax;
mcmcInfo.MS2_kernel = flipud(mcmcInfo.coeff_MS2(1:nStepsMax,:));

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
% index_helper1 = mcmcInfo.chain_id_ref*seq_len_temp + mcmcInfo.trace_id_ref*seq_len_temp*n_chains;
index_helper2 = mcmcInfo.chain_id_ref*nStates^2 + (mcmcInfo.row_ref-1)*nStates;
% index_helper3 = mcmcInfo.chain_id_ref*seq_len_dummy + mcmcInfo.trace_id_ref*seq_len_dummy*n_chains;

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
    
    %%% calculate fluorescence probability component
    % calculate fluo error term    

    logL_fluo = calculate_fluo_logL_v4(mcmcInfo);              
    
    %%% put everything together
    total_log_likelihoods = logL_tr + logL_fluo;    
    
    % apply differential temperature correction if appropirate
%     if mcmcInfo.temperingFlag 
%         total_log_likelihoods = total_log_likelihoods ./ mcmcInfo.tempGradVec;
%     end
    total_log_likelihoods = total_log_likelihoods - logsumexp(total_log_likelihoods,1);
    total_likelihoods = exp(total_log_likelihoods);    

    %%% draw new samples
    option_array = cumsum(total_likelihoods);
    rand_array = repmat(rand(1,n_chains,n_traces),nStates,1);
    sample_chains_temp(prev_time_index+1,:,:) = sum(rand_array > option_array) + 1;
    
    % pass along to dummy array    
    mcmcInfo.sample_chains_temp = sample_chains_temp;
   
end 

mcmcInfo.sample_chains = sample_chains_temp(2:end-1,:,:);
em_time = toc
if mcmcInfo.em_timer_flag
    mcmcInfo.em_time_vec(mcmcInfo.step) = em_time;
end
