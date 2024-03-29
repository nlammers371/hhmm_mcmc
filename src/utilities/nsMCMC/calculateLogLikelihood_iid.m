function mcmcInfo = calculateLogLikelihood_iid(mcmcInfo)
    
    pi0_log = log(mcmcInfo.pi0_curr);
    sigma_curr = mcmcInfo.sigma_curr;
%     nStates = size(pi0_curr,2);  
    n_traces = mcmcInfo.n_traces;
    n_chains = mcmcInfo.n_chains_eff;

    %%% record transitions

%     % get linear indices
%     from_array = mcmcInfo.sample_chains(1:end-1,:,:);
%     to_array = mcmcInfo.sample_chains(2:end,:,:);        
%     row_col_array = (from_array-1)*nStates+ to_array;% + nStates^2*mcmcInfo.chain_id_ref;
%     lin_index_array = row_col_array;
% 
%     % get state counts
%     for k = 1:nStates
%         mcmcInfo.initial_state_array(k,:) = sum(mcmcInfo.sample_chains(1,:,:)==k,3);% / n_chains;
%     end
% 
%     % get transition counts
%     unique_indices = unique(lin_index_array(:));    
% 
%     n_vec = (0:n_chains-1)*nStates^2;
%     for i = 1:length(unique_indices)        
%         mcmcInfo.transition_count_array(n_vec+unique_indices(i)) = sum(sum(lin_index_array==i,1),3);
%     end
% 
%     % calculate transition likelihoods
%     logL_transition_array = A_log(lin_index_array);

    % calculate fluo-based likelihood
    ref_fluo = repmat(permute(mcmcInfo.observed_fluo,[1 3 2]),1,n_chains,1);
    sigma_ref = repmat(sigma_curr',1,1,n_traces);    
    logL_fluo = -0.5*(((ref_fluo-mcmcInfo.sample_fluo)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));

    % calculate initial state likelihoods    
    logL_pi0 = pi0_log(mcmcInfo.sample_chains);    

    % record
    mcmcInfo.trace_logL_array = logL_pi0 + logL_fluo;
    mcmcInfo.trace_logL_vec = mean(mcmcInfo.trace_logL_array);

    update_flag = mod(mcmcInfo.step,mcmcInfo.update_increment) == 0 || mcmcInfo.step == 1 || mcmcInfo.step == mcmcInfo.n_mcmc_steps;
    update_index = floor((mcmcInfo.step-1)/mcmcInfo.update_increment) + 1;
    if mcmcInfo.step == 1
        update_index = 1;
    end
    if update_flag
        mcmcInfo.logL_vec(update_index,:) = mean(mcmcInfo.trace_logL_vec,3);%
    end