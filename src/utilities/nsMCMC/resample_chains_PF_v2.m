function mcmcInfo = resample_chains_PF_v2(mcmcInfo)

% This script resamples the microscopic promoter state for each extant
% chain using particle filtering methods that mimic classical
% forward-backward algorithm

% extract parameters
A = mcmcInfo.A_curr(:,:,1);
pi0 = mcmcInfo.pi0_curr(1,:);
nStates = mcmcInfo.nStates;
n_traces = mcmcInfo.n_traces;
n_chains = mcmcInfo.n_chains_eff;
nStepsMax = mcmcInfo.nStepsMax;
seq_len = mcmcInfo.seq_length;

% generate reference vector
mcmcInfo.chain_id_ref = 0:n_chains-1;
mcmcInfo.trace_id_ref = reshape(0:n_traces-1,1,1,[]);
mcmcInfo.row_ref = (1:nStates)';
mcmcInfo.step_ref = (-nStepsMax+1:nStepsMax-1)';

% generate temporary chain array that includes post and prior states to
% make resampling easier below
sample_chains_fwd = NaN(size(mcmcInfo.sample_chains));
% sample_emissions_temp = NaN(size(mcmcInfo.sample_chains));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conduct forward particle simulations first
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial_state_vec = NaN(1,n_chains,n_traces);
% for n = 1:n_chains
%     initial_state_vec(1,n,:) = reshape(randsample(1:nStates,n_traces,true,pi0(n,:)),1,1,n_traces);    
% end 

% generate fluorescent noise lookup table 
sigma_ref = repmat(mcmcInfo.sigma_curr(1),1,n_chains,n_traces);
fluo_ref = permute(mcmcInfo.observed_fluo,[1,3,2]);

% generate lookup table for fluorescence values corresponding to each state
% at each position along the gene body
max_w = find(mcmcInfo.coeff_MS2(:,1)>0,1,'last');
ms2_kernel = mcmcInfo.coeff_MS2(1:max_w,:);

% generate array for conversion from longform to short CP state notation
cp_vec = (nStates).^(max_w-1:-1:0);
cp_conv_array = repmat(cp_vec',1,n_chains,n_traces);

% initialize array to keep track of compound states
cp_state_tracker = ones(max_w,n_chains,n_traces);
cp_init_tracker = zeros(max_w,n_chains,n_traces);
cp_state_array_fwd = zeros(seq_len,n_chains,n_traces);        
F_array_fwd_long = zeros(seq_len,n_chains,n_traces,nStates);        

% array to store cp state fluo likelihoods
sample_fluo_prob_fwd = zeros(seq_len,n_chains,n_traces);        

% initialize new "predicted fluo" array
sample_fluo_fwd = NaN(size(sample_chains_fwd));
sample_logL_fwd = NaN(size(sample_chains_fwd));

% initialize tracker for state transitions
% tr_tracker_fwd = zeros(nStates, nStates, seq_len-1, n_traces);

tic

% iterate through indices to sample
for i = 1:seq_len
    
    mcmcInfo.samp_index = i;        
    prev_time_index = i - 1;
    
    %%% simulate next state transition %%%
    
    % calculate linear indices 
    if i > 1
        prev_state_array = sample_chains_fwd(prev_time_index,:,:);
    
        row_col_array_from = ...
                             (prev_state_array-1)*nStates + mcmcInfo.row_ref;    % mcmcInfo.chain_id_ref*nStates^2 + ...
    
        % extract probabilities
        tr_probs = A(row_col_array_from);
    
        % simulate jumps
        option_array = cumsum(tr_probs);
        rand_array_prop = repmat(rand(1,n_chains,n_traces),nStates,1);
        new_states = sum(rand_array_prop > option_array) + 1;
                        
        % calculate likelihood associated with each choice
        new_state_tr_vec = log(tr_probs(new_states+nStates*mcmcInfo.chain_id_ref));
        
    else                
        new_states = reshape(randsample(1:nStates,n_traces*n_chains,true,pi0),1,n_chains,n_traces);            
    end    
    sample_chains_fwd(i,:,:) = new_states;  
    
    
    
    %%% update compound state tracker
    cp_state_tracker = circshift(cp_state_tracker,1,1);
    cp_state_tracker(1,:,:) = new_states;        
       
    %%% initiation state tracker
    cp_init_tracker = circshift(cp_init_tracker,1,1);
    v_lin_indices = (mcmcInfo.chain_id_ref+1) + n_chains*(new_states-1);
    cp_init_tracker(1,:,:) = mcmcInfo.v_curr(v_lin_indices);
    
    %%% calculate fluorescence
    fluo_next = sum(cp_init_tracker.*ms2_kernel,1);    
    
    %%% resample based on difference from observed fluo %%%
    
    % calculate fluo probabilities 
    logL_fluo_full = -0.5*(((reshape(fluo_ref(i,:),1,1,[])-fluo_next)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));
    probs_fluo = exp(logL_fluo_full-logsumexp(logL_fluo_full,2));    
       
    % resample
    option_array = [zeros(1,1,n_traces) cumsum(probs_fluo,2)];
    rand_array = rand(1,n_chains,n_traces);    
    
    for t = 1:n_traces
      
        % resample
        rs_ids = discretize(rand_array(1,:,t),option_array(1,:,t));
        
        % resample arrays
        sample_chains_fwd(i,:,t) = sample_chains_fwd(i,rs_ids,t);
        sample_fluo_fwd(i,:,t) = fluo_next(1,rs_ids,t);
        if i > 1
            prev_state_array(1,:,t) = prev_state_array(1,rs_ids,t);
        end
        
        % resample compound state trackers
        cp_state_tracker(:,:,t) = cp_state_tracker(:,rs_ids,t);
        cp_init_tracker(:,:,t) = cp_init_tracker(:,rs_ids,t);
        
        % add fluo probabilities
        sample_fluo_prob_fwd(i,:,t) = probs_fluo(1, rs_ids, t);
    end
    
    % convert cp states to integer representation   
    cp_state_array_fwd(i,:,:) = 1 + sum(cp_conv_array.*(cp_state_tracker-1),1);    
    for k = 1:nStates
%         test(k,:,:) = sum(ms2_kernel.*(cp_state_tracker==k),1);
        F_array_fwd_long(i,:,:,k) = sum(ms2_kernel.*(cp_state_tracker==k),1);
    end
%     if i > 1
%         tr_counts_temp = tr_tracker_fwd(:,:,i-1,:);
%         for n = 1:nStates
%             for m = 1:nStates
%                 tr_counts_temp(n,m,1,:) = reshape(sum(prev_state_array==m & sample_chains_fwd(i,:,t)==n,2),1,1,1,[]);
%             end
%         end
%         tr_tracker_fwd(:,:,i-1,:) = tr_counts_temp;
%     end
end 
% toc


%%
em_time = toc;
if mcmcInfo.em_timer_flag
    mcmcInfo.em_time_vec(mcmcInfo.step) = em_time;
end