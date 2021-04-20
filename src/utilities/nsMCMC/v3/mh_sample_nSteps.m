function mcmcInfo = mh_sample_nSteps(mcmcInfo)

% pull out parameters for conveneice
n_chains = mcmcInfo.n_chains;

% iterate through traces
for m = 1:n_traces
  
    % extract relevant slice
    chain_slice = mcmcInfo.sample_chains(:,:,m);

    % extract trace
    ref_trace = mcmcInfo.observed_fluo(:,m);          

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Propose new W values 
    %%%%%%%%%%%%%%%%%%%%%%%%%                               
    linear_indices = (trace_array_full-1)*n_chains_eff + chain_id_array_full+1;
    initiation_rates = mcmcInfo.v_curr(linear_indices);    
    
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