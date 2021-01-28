function mcmcInfo = get_empirical_counts_ens(mcmcInfo)
    
    % extract parameters
    A_curr = mcmcInfo.A_curr;
    A_log = log(mcmcInfo.A_curr);
    pi0_log = log(mcmcInfo.pi0_curr);
    sigma_curr = mcmcInfo.sigma_curr;
    nStates = size(A_curr,1);
    nSteps = mcmcInfo.nSteps;
    n_traces = mcmcInfo.n_traces;
    n_chains = mcmcInfo.n_chains;
    seq_length = mcmcInfo.seq_length;
    
    % generate reference vector
%     chain_id_ref = 0:n_chains-1;

    % initialize count arrays 
    mcmcInfo.transition_count_array = zeros(size(A_curr,1),size(A_curr,2));        
    mcmcInfo.initial_state_array = zeros(nStates,1);   
    
    % initialize likelihood array
    mcmcInfo.trace_logL_array = -Inf(seq_length,n_chains,n_traces);
    mcmcInfo.trace_logL_vec = -Inf(1,n_chains,n_traces);
    
        
    %%% record transitions
    
    % get linear indices
    from_array = mcmcInfo.sample_chains(1:end-1,:,:);
    to_array = mcmcInfo.sample_chains(2:end,:,:);        
    row_col_array = (from_array-1)*nStates+to_array;
    lin_index_array = row_col_array;
    
    % get state counts
    for k = 1:nStates
        mcmcInfo.initial_state_array(k) = sum(sum(mcmcInfo.sample_chains(1,:,:)==k,3)) / n_chains;
    end
    
    % get transition counts
    unique_indices = unique(lin_index_array(:));
    transition_counts = histc(lin_index_array(:),unique_indices);
    
    for i = 1:length(unique_indices)
        mcmcInfo.transition_count_array(unique_indices(i)) = transition_counts(i) / n_chains;
    end
    
    % calculate transition likelihoods
    logL_transition_array = A_log(lin_index_array);
    
    % calculate fluo-based likelihood
    ref_fluo = repmat(permute(mcmcInfo.observed_fluo,[1 3 2]),1,n_chains,1);
    sigma_ref = repmat(sigma_curr',size(ref_fluo,1),1,n_traces);    
    logL_fluo = -0.5*(((ref_fluo-mcmcInfo.sample_fluo)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));
        
    % calculate initial state likelihoods    
    logL_pi0 = pi0_log(mcmcInfo.sample_chains(1,:,:));  
    if n_chains == 1
        logL_pi0 = reshape(logL_pi0,1,1,[]);
    end
    logL_transition_array_full = cat(1,logL_pi0,logL_transition_array);
    
    % record
    mcmcInfo.trace_logL_array = logL_transition_array_full + logL_fluo;
    mcmcInfo.trace_logL_vec = mean(mcmcInfo.trace_logL_array);
    
    update_flag = mod(mcmcInfo.step-1,mcmcInfo.update_increment) == 0 || mcmcInfo.step-1 == 1;
    update_index = (mcmcInfo.step-1)/mcmcInfo.update_increment + 1;
    if mcmcInfo.step-1 == 1
        update_index = 1;
    end
    if update_flag
        mcmcInfo.logL_vec(update_index,:) = mean(mean(mcmcInfo.trace_logL_vec,3));%
    end
    