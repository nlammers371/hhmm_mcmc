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

% generate block to store state guesses
% sim_block = NaN(2*us_factor+1,n_chains,n_traces);

% generate helper arrays with pre-calculated transition probs
% tr_prob_cell = cell(nStates,nStates,n_chains);
% iter_vec_fw = reshape(1:us_factor,1,1,[]);
% iter_vec_bk = reshape(us_factor:-1:1,1,1,[]);


% for n = 1:n_chains
%     A_temp = A(:,:,n);
%     A_array_fw = NaN(nStates,nStates,us_factor);
%     A_array_bk = NaN(nStates,nStates,us_factor);
%     for u = 1:us_factor
%         A_array_fw(:,:,u) = A_temp^iter_vec_fw(u);
%         A_array_bk(:,:,u) = A_temp^iter_vec_bk(u);
%     end
%     for k = 1:nStates
%         for j = 1:nStates
%             tr_slice = permute(A_array_fw(:,j,:),[1 3 2]).*permute(A_array_bk(k,:,:),[2 3 1]);
%             tr_slice = tr_slice ./ sum(tr_slice,1);
%             tr_prob_cell{k,j,n} = tr_slice;
%         end
%     end                      
% end  

% sample_indices = 
% iterate through indices to sample
for i = [(1:seq_len) ((seq_len-1):-1:1)]
    
    mcmcInfo.samp_index = i;%sample_indices(i);        
    prev_time_index = i;%sample_indices(i);    
    post_time_index = prev_time_index + 2;
    
    %%% previous state %%%
%     samp_index_micro = us_factor*mcmcInfo.samp_index;
    prev_time_index_micro = (i-1)*us_factor+1;
    post_time_index_micro = i*us_factor+2;
    try
      sim_block = mcmcInfo.sample_chains_temp(prev_time_index_micro:post_time_index_micro,:,:);
    catch
      error('afa')
    end
    mcmcInfo.relevant_indices = prev_time_index_micro+1:post_time_index_micro-1;
    mcmcInfo.curr_block = sim_block(2:end-1,:,:);
%     mcmcInfo.fluo_curr_block = sample_fluo_temp(prev_time_index_micro-us_factor+2:post_time_index_micro-us_factor,:,:);
    
    % conduct forward simulation    
    for t = [2:us_factor-1 us_factor-1:-1:2]
      
        % calculate linear indices 
        prev_state_array = sim_block(t-1,:,:);
        row_col_array_from = mcmcInfo.chain_id_ref*nStates^2 + (prev_state_array-1)*nStates+mcmcInfo.row_ref;    
    
        % extract probabilities
        prev_probs = A(row_col_array_from);
          
        % calculate linear indices 
        post_state_array = sim_block(t+1,:,:);
        row_col_array_to = index_helper2 + post_state_array;
        
        % extract probabilities                
        post_probs = A(row_col_array_to);
        post_probs = post_probs ./ sum(post_probs,1);
        
        tr_probs = post_probs.*prev_probs;
        tr_probs = tr_probs ./ sum(tr_probs);
        
        %%% draw new samples
        option_array = cumsum(tr_probs);
        rand_array = repmat(rand(1,n_chains,n_traces),nStates,1);
        sim_block(t,:,:) = sum(rand_array > option_array) + 1;
    end
                   
    mcmcInfo.prop_block = sim_block(2:end-1,:,:);    
    
    %%% calculate fluorescence probability component
    % calculate fluo error term        
    mcmcInfo = calculate_fluo_logL_v5_MH(mcmcInfo);        
   
end 

mcmcInfo.sample_chains = mcmcInfo.sample_chains_temp(2:end-1,:,:);
mcmcInfo.sample_fluo = mcmcInfo.sample_fluo_temp;
mcmcInfo = rmfield(mcmcInfo,{'sample_chains_temp','sample_chains_temp','sample_chains_temp'});

if mcmcInfo.em_timer_flag
    em_time = toc;
    mcmcInfo.em_time_vec(mcmcInfo.step) = em_time;
end
