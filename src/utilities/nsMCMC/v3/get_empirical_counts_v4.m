function mcmcInfo = get_empirical_counts_v4(mcmcInfo)
    
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
    mcmcInfo.transition_count_array = zeros(size(A_curr,1),size(A_curr,2),n_chains);    
    mcmcInfo.initial_state_array = zeros(nStates,n_chains);   
    
    % initialize likelihood array
%     mcmcInfo.trace_logL_array = -Inf(seq_length,n_chains,n_traces);
    mcmcInfo.trace_logL_vec = NaN(n_traces,n_chains);
    mcmcInfo.trace_n_points = NaN(1,n_traces);
            
    %%% record transitions
    
    % get linear indices
    for m = 1:n_traces
        from_array = mcmcInfo.sample_chains{m}(1:end-1,:,1);
        to_array = mcmcInfo.sample_chains{m}(2:end,:,1);        
        row_col_array = (from_array-1)*nStates+ to_array;% + nStates^2*mcmcInfo.chain_id_ref;
        lin_index_array = row_col_array;

        % get initial state counts
        for k = 1:nStates
            mcmcInfo.initial_state_array(k,:) = mcmcInfo.initial_state_array(k,:) + sum(mcmcInfo.sample_chains{m}(1,:,1)==k,3);% / n_chains;
        end

        % get transition counts
        unique_indices = unique(lin_index_array(:));    

        n_vec = (0:n_chains-1)*nStates^2;
        for i = 1:length(unique_indices)   
            ind_temp = n_vec+unique_indices(i);
            mcmcInfo.transition_count_array(ind_temp) = mcmcInfo.transition_count_array(ind_temp) + sum(sum(lin_index_array==i,1),3);
        end

        % calculate transition likelihoods
        logL_transition_array = A_log(lin_index_array);

        % calculate fluo-based likelihood
        ref_fluo = mcmcInfo.observed_fluo{m};
        sigma_ref = sigma_curr;%repmat(sigma_curr,size(ref_fluo,1),1,n_traces);    
        logL_fluo = -0.5*(((ref_fluo-mcmcInfo.sample_fluo{m}(:,:,1))./sigma_ref).^2 + log(2*pi*sigma_ref.^2));

        % calculate initial state likelihoods    
        logL_pi0 = pi0_log(mcmcInfo.sample_chains{m}(1,:,1));  
%         if n_chains == 1
%             logL_pi0 = reshape(logL_pi0,1,1,[]);
%         end
        logL_transition_array_full = cat(1,logL_pi0,logL_transition_array);

        % record
        trace_logL_array = logL_transition_array_full + logL_fluo;
        mcmcInfo.trace_logL_vec(m,:) = sum(trace_logL_array);
        mcmcInfo.trace_n_points(m) = size(trace_logL_array,1);
    end
    
    update_flag = mod(mcmcInfo.step-1,mcmcInfo.update_increment) == 0 || mcmcInfo.step-1 == 1 || mcmcInfo.step == mcmcInfo.n_mcmc_steps;
    update_index = floor((mcmcInfo.step-1)/mcmcInfo.update_increment) + 1;
    if mcmcInfo.step-1 == 1
        update_index = 1;
    end
    if update_flag
        logL_total = sum(mcmcInfo.trace_logL_vec)./sum(mcmcInfo.trace_n_points);
        mcmcInfo.logL_vec(update_index,:) = logL_total;%
    end
    