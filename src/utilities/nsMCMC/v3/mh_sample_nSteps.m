function mh_sample_nSteps(mcmcInfo)


% iterate through traces
for m = 1:n_traces
  
    % extract relevant slice
    chain_slice = mcmcInfo.sample_chains(:,:,m);

    % extract trace
    ref_trace = mcmcInfo.observed_fluo(:,m);          

%     %%%%%%%%%%%%%%%%%%%%%%%%%
%     % Transition component 
% 
%     %%% previous state %%%
%     % calculate linear indices 
%     prev_state_array = trace_array_full(swap_ind_array-1);
%     curr_state_array = trace_array_full(swap_ind_array);
%     row_col_array_from = chain_id_array_full.*nStates^2 + (prev_state_array-1)*nStates + curr_state_array;    
% 
%     % extract probabilities
%     prev_probs_log = A_log(row_col_array_from);
% 
%     %%% following state %%%
%     % calculate linear indices 
%     post_state_array = trace_array_full(swap_ind_array+1);
%     row_col_array_to = chain_id_array_full.*nStates^2 + (curr_state_array-1)*nStates + post_state_array;
% 
%     % extract probabilities
%     post_probs_log = A_log(row_col_array_to);
% 
%     % combine
%     logL_tr = prev_probs_log + post_probs_log;

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Emission component 
    %%%%%%%%%%%%%%%%%%%%%%%%%                               

    fluo_predicted = convn(coeff_MS2,initiation_rates,'full');             
    fluo_predicted = fluo_predicted(1:end-length(coeff_MS2)+1,:,:);

    sigma_ref = mcmcInfo.sigma_curr(chain_id_array_full+1);    

    logL_fluo = -0.5*(((ref_trace-fluo_predicted)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));


    % combine
    total_log_likelihoods = logL_tr + sum(logL_fluo);  

    % calculate MH metric
    L_factor = exp((total_log_likelihoods(1,:,3)-total_log_likelihoods(1,:,1))./temperature_array(1,:) + ...
                    (total_log_likelihoods(1,:,4)-total_log_likelihoods(1,:,2))./temperature_array(2,:));

    % perform MH move
    mh_flags = L_factor >= rand(size(L_factor));

    % update
    trace_array(swap_lin_indices(mh_flags)) = trace_array_prop(swap_lin_indices(mh_flags));
end