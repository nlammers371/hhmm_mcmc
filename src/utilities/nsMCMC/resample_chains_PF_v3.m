function mcmcInfo = resample_chains_PF_v3(mcmcInfo)

% This script resamples the microscopic promoter state for each extant
% chain using particle filtering methods that mimic classical
% forward-backward algorithm

% extract parameters
A = mcmcInfo.A_curr;
pi0 = mcmcInfo.pi0_curr;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conduct forward particle simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize array to store forward sweep results
sample_chains_fwd = NaN(size(mcmcInfo.sample_chains));

% initialize new "predicted fluo" array
sample_fluo_fwd = NaN(size(sample_chains_fwd));

% generate fluorescent noise lookup table 
sigma_ref = repmat(mcmcInfo.sigma_curr',1,1,n_traces);
fluo_ref = permute(mcmcInfo.observed_fluo,[1,3,2]);

% generate lookup table with ms2 kernel values
max_w = 0;
for n = 1:n_chains
    max_w = max([max_w find(mcmcInfo.coeff_MS2(:,n)>0,1,'last')]);
end    
ms2_kernel = mcmcInfo.coeff_MS2(1:max_w,:);

% generate array for conversion from longform to short CP state notation
% cp_vec = (nStates).^(max_w-1:-1:0);
% cp_conv_array = repmat(cp_vec',1,n_chains,n_traces);

% initialize array to keep track of compound states
cp_state_tracker = ones(max_w,n_chains,n_traces);
cp_init_tracker = zeros(max_w,n_chains,n_traces);
% cp_state_array_fwd = zeros(seq_len,n_chains,n_traces);                 

% tic
% iterate through indices to sample
for i = 1:seq_len
    
    mcmcInfo.samp_index = i;        
    prev_time_index = i - 1;
    
    %%% simulate next state transition %%%
    
    % calculate linear indices 
    if i > 1
        prev_state_array = sample_chains_fwd(prev_time_index,:,:);
    
        row_col_array_from = mcmcInfo.chain_id_ref.*nStates^2 + ...
                             (prev_state_array-1)*nStates + mcmcInfo.row_ref;    % mcmcInfo.chain_id_ref*nStates^2 + ...
    
        % extract probabilities
        tr_probs = log(A(row_col_array_from));
                        
    else                
        tr_probs = log(repmat(pi0',1,1,n_traces));%new_states = reshape(randsample(1:nStates,n_traces*n_chains,true,pi0),1,n_chains,n_traces);            
    end
    
    %%% calculate predicted fluo values corresponding to each possible
    %%% transition    
    fluo_probs = NaN(size(tr_probs));
    cp_init_tracker = circshift(cp_init_tracker,1,1);
    for k = 1:nStates
        prop_states = k*ones(1,n_chains,n_traces);
        v_lin_indices = (mcmcInfo.chain_id_ref+1) + n_chains*(prop_states-1);
        cp_init_tracker(1,:,:) = mcmcInfo.v_curr(v_lin_indices);
        
        % calculate fluorescence
        fluo_prop = sum(cp_init_tracker.*ms2_kernel,1); 
        
        % calculate fluo probabilities 
        fluo_probs(k,:,:) = -0.5*(((fluo_ref(i,:,:)-fluo_prop)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));
    end
    
%     fluo_probs = fluo_probs - logsumexp(fluo_probs,1);
    
    %%% combine probabilities and sample
    log_probs_total = fluo_probs + tr_probs;
    probs_total = exp(log_probs_total - logsumexp(log_probs_total,1));
    
    option_array = cumsum(probs_total,1);
    rand_array_prop = repmat(rand(1,n_chains,n_traces),nStates,1);
    new_states = sum(rand_array_prop > option_array) + 1;

    %%% update compound state tracker
    cp_state_tracker = circshift(cp_state_tracker,1,1);
    cp_state_tracker(1,:,:) = new_states;                 
    
    %%% update initiation state tracker
    v_lin_indices = (mcmcInfo.chain_id_ref+1) + n_chains*(new_states-1);
    cp_init_tracker(1,:,:) = mcmcInfo.v_curr(v_lin_indices);
    
    %%% update fluo and promoter state arrays
    sample_chains_fwd(i,:,:) = new_states;
    sample_fluo_fwd(i,:,:) = sum(cp_init_tracker.*ms2_kernel,1); 
end 
% toc

mcmcInfo.sample_chains = sample_chains_fwd;
mcmcInfo.sample_fluo = sample_fluo_fwd;

%%
em_time = toc;
if mcmcInfo.em_timer_flag
    mcmcInfo.em_time_vec(mcmcInfo.step) = em_time;
end