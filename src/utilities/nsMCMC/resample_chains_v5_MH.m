function mcmcInfo = resample_chains_v5_MH(mcmcInfo)

if mcmcInfo.em_timer_flag
    tic
end
% tic

% extract parameters
A = mcmcInfo.A_curr;
pi0 = mcmcInfo.pi0_curr;
nStates = mcmcInfo.nStates;
n_traces = mcmcInfo.n_traces;
n_chains = mcmcInfo.n_chains_eff;
us_factor = mcmcInfo.upsample_factor;

% nStepsMax = mcmcInfo.nStepsMax;
seq_len = mcmcInfo.seq_length;%mcmcInfo.seq_length;

% generate lookup table with ms2 kernel values
nStepsMax = 0;
for n = 1:n_chains
    nStepsMax = max([nStepsMax find(mcmcInfo.coeff_MS2(:,n)>0,1,'last')]);
end    
mcmcInfo.max_w = nStepsMax;
mcmcInfo.MS2_kernel = mcmcInfo.coeff_MS2(1:nStepsMax,:);
mcmcInfo.MS2_kernel_us = mcmcInfo.coeff_MS2_us(1:nStepsMax*us_factor,:);

% note that this version uses normal array indexing
mcmcInfo.v_curr_long = repmat(permute(mcmcInfo.v_curr',[3 2 1]),nStepsMax*us_factor,1,1);

% generate reference vector
mcmcInfo.chain_id_ref = 0:n_chains-1;
mcmcInfo.trace_id_ref = reshape(0:n_traces-1,1,1,[]);
mcmcInfo.row_ref = (1:nStates)';

% generate random sampling orders
% n_reps = mcmcInfo.n_reps;
% sample_indices = randsample(repelem(1:seq_len,n_reps),n_reps*seq_len,false);

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
sample_fluo_temp = mcmcInfo.sample_fluo;
mcmcInfo.sample_fluo_temp = sample_fluo_temp;

% preallocate helper array for indexing
index_helper2 = mcmcInfo.chain_id_ref*nStates^2 + (mcmcInfo.row_ref-1)*nStates;

sim_block = NaN(2*us_factor+1,n_chains,n_traces);

% sample_indices = 
% iterate through indices to sample
for i = [(1:seq_len) ((seq_len-1):-1:1)]
    
    mcmcInfo.samp_index = i;%sample_indices(i);        
    prev_time_index = i;%sample_indices(i);    
    post_time_index = prev_time_index + 2;
    
    %%% previous state %%%
    samp_index_micro = us_factor*mcmcInfo.samp_index;
    prev_time_index_micro = us_factor*prev_time_index;
    post_time_index_micro = us_factor*post_time_index;
    
    sim_block(1,:,:) = sample_chains_temp(prev_time_index_micro,:,:);
    sim_block(end,:,:) = sample_chains_temp(post_time_index_micro,:,:);
    
    % conduct forward simulation    
    for t = 1:us_factor-1
      
        % calculate linear indices 
        prev_state_array = sim_block(t,:,:);
        row_col_array_from = mcmcInfo.chain_id_ref*nStates^2 + (prev_state_array-1)*nStates+mcmcInfo.row_ref;    
    
        % extract probabilities
        prev_probs = A(row_col_array_from);
          
        %%% draw new samples
        option_array = cumsum(prev_probs);
        rand_array = repmat(rand(1,n_chains,n_traces),nStates,1);
        sim_block(t+1,:,:) = sum(rand_array > option_array) + 1;
    end
        
    % conduct backward simulation      
    for t = 2*us_factor+1:-1:us_factor+1
      
        % calculate linear indices 
        post_state_array = sim_block(t,:,:);
        row_col_array_to = index_helper2 + post_state_array;
        
        % extract probabilities                
        post_probs = A(row_col_array_to);
        post_probs = post_probs ./ sum(post_probs,1);
                     
        %%% draw new samples
        option_array = cumsum(post_probs);
        rand_array = repmat(rand(1,n_chains,n_traces),nStates,1);
        sim_block(t-1,:,:) = sum(rand_array > option_array) + 1;
    end
    
    % now draw the middle state
    % calculate linear indices 
    prev_state_array = sample_chains_temp(samp_index_micro-1,:,:);
    row_col_array_from = mcmcInfo.chain_id_ref*nStates^2 + (prev_state_array-1)*nStates+mcmcInfo.row_ref;    
    
    % extract probabilities
    prev_probs = A(row_col_array_from);
          
    %%% following state %%%
    % calculate linear indices 
    post_state_array = sample_chains_temp(samp_index_micro+1,:,:);  
    row_col_array_to = index_helper2 + post_state_array;
        
    % extract probabilities
    post_probs = A(row_col_array_to);        
       
    % combine
    prob_tr = prev_probs.*post_probs;
    %%% draw new samples
    option_array = cumsum(prob_tr);
    rand_array = repmat(rand(1,n_chains,n_traces),nStates,1);
    sim_block(us_factor+1,:,:) = sum(rand_array > option_array) + 1;
    mcmcInfo.sim_block = sim_block;
    
    %%% calculate fluorescence probability component
    % calculate fluo error term        
    [logL_fluo, new_fluo_array, comp_indices] = calculate_fluo_logL_v5_MH(mcmcInfo);       
 
    %%% put everything together
    total_log_likelihoods = logL_tr + logL_fluo;    
    
    total_log_likelihoods = total_log_likelihoods - logsumexp(total_log_likelihoods,1);
    total_likelihoods = exp(total_log_likelihoods);    

    %%% draw new samples
    option_array = cumsum(total_likelihoods);
    rand_array = repmat(rand(1,n_chains,n_traces),nStates,1);
    new_states = sum(rand_array > option_array) + 1;
    sample_chains_temp(prev_time_index+1,:,:) = new_states;
    
    % update fluo_array
    fluo_update_indices_to = comp_indices + seq_len*mcmcInfo.chain_id_ref + n_chains*seq_len*mcmcInfo.trace_id_ref;
    n_indices = length(comp_indices);
    fluo_update_indices_from = (1:n_indices)' + n_indices*mcmcInfo.chain_id_ref + n_chains*n_indices*mcmcInfo.trace_id_ref + n_chains*n_indices*n_traces*(new_states-1);    
    mcmcInfo.sample_fluo_temp(fluo_update_indices_to) = new_fluo_array(fluo_update_indices_from);
    
    % pass along to dummy array    
    mcmcInfo.sample_chains_temp = sample_chains_temp;
   
end 

mcmcInfo.sample_chains = sample_chains_temp(2:end-1,:,:);
mcmcInfo.sample_fluo = mcmcInfo.sample_fluo_temp;
mcmcInfo = rmfield(mcmcInfo,{'sample_chains_temp','sample_chains_temp','sample_chains_temp'});

if mcmcInfo.em_timer_flag
    em_time = toc;
    mcmcInfo.em_time_vec(mcmcInfo.step) = em_time;
end
