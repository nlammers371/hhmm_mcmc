function mcmcInfo = resample_chains_PF_LL(mcmcInfo)

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

% generate fluorescent noise lookup table 
sigma_ref = repmat(mcmcInfo.sigma_curr(1),1,n_chains,n_traces);
fluo_ref = permute(mcmcInfo.observed_fluo,[1,3,2]);

% generate lookup table for fluorescence values corresponding to each state
% at each position along the gene body
max_w = find(mcmcInfo.coeff_MS2(:,1)>0,1,'last');
ms2_kernel = mcmcInfo.coeff_MS2(1:max_w,:);

% generate array for conversion from longform to short CP state notation
% cp_vec = (nStates).^(max_w-1:-1:0);
% cp_conv_array = repmat(cp_vec',1,n_chains,n_traces);

% initialize array to keep track of compound states
cp_state_tracker = ones(max_w,n_chains,n_traces);
cp_init_tracker = zeros(max_w,n_chains,n_traces);
% cp_state_array_fwd = zeros(seq_len,n_chains,n_traces);        
% F_array_fwd_long = zeros(seq_len,n_chains,n_traces,nStates);        

% array to store cp state fluo likelihoods
% sample_fluo_prob_fwd = zeros(seq_len,n_chains,n_traces);        

% initialize new "predicted fluo" array
sample_fluo_fwd = NaN(size(sample_chains_fwd));
sample_logL_fwd = NaN(size(sample_chains_fwd));
sample_parent_fwd = NaN(size(sample_chains_fwd,1)+1,n_chains,n_traces);
sample_parent_fwd(1,:,:) = repmat(1:n_chains,1,1,n_traces);
sample_parent_cum_fwd = NaN(size(sample_chains_fwd,1)+1,n_chains,n_traces);
sample_parent_cum_fwd(1,:,:) = repmat(1:n_chains,1,1,n_traces);

% initialize tracker for state transitions
% tr_tracker_fwd = zeros(nStates, nStates, seq_len-1, n_traces);

% tic

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
        new_state_tr_vec = log(pi0(new_states));
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
    sample_fluo_fwd(i,:,:) = fluo_next;
    %%% resample based on difference from observed fluo %%%
    
    % calculate chain probabilities 
    logL_fluo_full = -0.5*(((reshape(fluo_ref(i,:),1,1,[])-fluo_next)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));    
    
%         fluo_logL_fwd(i,:,:) = fluo_logL_fwd(i-1,:,:) + logL_fluo_full;
    sample_logL_fwd(i,:,:) = logL_fluo_full + new_state_tr_vec;
%         fluo_logL_fwd(i,:,:) = logL_fluo_full;
    
    % calculate RS factor
    rs_factor = mean(exp(2*logsumexp(sample_logL_fwd(i,:,:),2) - logsumexp(2*sample_logL_fwd(i,:,:),2)),3);
    
    if rs_factor < n_chains*0.5
      
        % check to see if we need to resample
        probs_chain = exp(sample_logL_fwd(i,:,:)-logsumexp(sample_logL_fwd(i,:,:),2));    

        % resample
        option_array = [zeros(1,1,n_traces) cumsum(probs_chain,2)];
        rand_array = rand(1,n_chains,n_traces);    
    
        for t = 1:n_traces

            % resample
            rs_ids = discretize(rand_array(1,:,t),option_array(1,:,t));

            % resample arrays
            sample_chains_fwd(i,:,t) = sample_chains_fwd(i,rs_ids,t);
            sample_fluo_fwd(i,:,t) = sample_fluo_fwd(i,rs_ids,t);
            sample_logL_fwd(i,:,t) = sample_logL_fwd(i,rs_ids,t);
            sample_parent_fwd(i+1,:,t) = rs_ids;
            sample_parent_cum_fwd(i+1,:,t) = sample_parent_cum_fwd(i,rs_ids,t);
            if i > 1
                prev_state_array(1,:,t) = prev_state_array(1,rs_ids,t);
            end

            % resample compound state trackers
            cp_state_tracker(:,:,t) = cp_state_tracker(:,rs_ids,t);
            cp_init_tracker(:,:,t) = cp_init_tracker(:,rs_ids,t);

        end
%         fluo_logL_fwd(i,:,:) = 0;
    else
        sample_parent_fwd(i+1,:,:) = repmat(1:n_chains,1,1,n_traces);
        sample_parent_cum_fwd(i+1,:,:) = sample_parent_cum_fwd(i,:,:);
    end
   
end

mcmcInfo.sample_fluo = sample_fluo_fwd;
mcmcInfo.sample_chains = sample_chains_fwd;
mcmcInfo.sample_logL = sample_logL_fwd;
mcmcInfo.LL = 5*mean(mean(logsumexp(sample_logL_fwd,2)-log(n_chains),3),1);% - log(n_chains*n_traces);
% toc

