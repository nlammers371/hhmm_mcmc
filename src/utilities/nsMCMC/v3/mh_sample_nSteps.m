function mcmcInfo = mh_sample_nSteps(mcmcInfo)

% pull out parameters for conveneice
n_chains = mcmcInfo.n_chains;
n_traces = mcmcInfo.n_traces;
n_chains_eff = mcmcInfo.n_chains_eff;
chain_id_eff_ref = 0:n_chains_eff-1;

% initialize array to store log likelihoods 
trace_logL_array = NaN(n_traces,n_chains_eff);
trace_temp_array = NaN(size(mcmcInfo.sample_fluo));

% propose new nSteps values
nStepsCurr = mcmcInfo.nStepsCurr;
ubVec = mcmcInfo.nStepsMax-nStepsCurr;
lbVec = mcmcInfo.nStepsMin-nStepsCurr;
nStepsProp = nStepsCurr + trandn(lbVec,ubVec)*mcmcInfo.nStepsPropSize;

% generate MS2 kernels
coeff_MS2_prop = NaN(size(mcmcInfo.coeff_MS2));
for n = 1:length(nStepsProp)
    alpha = mcmcInfo.alpha_frac*nStepsProp(n);
    coeff_MS2_prop(:,n) = ms2_loading_coeff_frac(alpha, nStepsProp(n), mcmcInfo.nStepsMax);
end

% iterate through traces
for m = 1:n_traces
  
    % extract relevant slice
    chain_slice = mcmcInfo.sample_chains(:,:,m);

    % extract trace
    ref_trace = mcmcInfo.observed_fluo(:,m);          

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Propose new W values 
    %%%%%%%%%%%%%%%%%%%%%%%%%                               
    linear_indices = (chain_slice-1)*n_chains_eff + chain_id_eff_ref + 1;
    initiation_rates = mcmcInfo.v_curr(linear_indices);    
    if n_chains_eff==1
        initiation_rates = initiation_rates';
    end
    % extract current predicted fluo
    fluo_curr = mcmcInfo.sample_fluo(:,:,m);        
    
    % get predicted fluo    
    fluo_prop = fluo_conv_fun(initiation_rates,coeff_MS2_prop);
    sigma_ref = mcmcInfo.sigma_curr';    

    logL_fluo = -0.5*(((ref_trace-cat(3,fluo_prop,fluo_curr))./sigma_ref).^2 + log(2*pi*sigma_ref.^2));
    
    % reshape and store difference
    total_log_likelihoods = permute(sum(logL_fluo),[3 2 1]);      
    trace_logL_array(m,:) = (total_log_likelihoods(1,:)-total_log_likelihoods(2,:))./mcmcInfo.tempGradVec;
    
    % store fluo
    trace_temp_array(:,:,m) = fluo_prop;
end

% calculate MH metric
L_factor = exp(sum(trace_logL_array));

% perform MH move
mh_flags = L_factor >= rand(size(L_factor));

% update    
mcmcInfo.coeff_MS2(:,mh_flags) = coeff_MS2_prop(:,mh_flags);
mcmcInfo.nStepsCurr(mh_flags) = nStepsProp(mh_flags);
mcmcInfo.sample_fluo(:,mh_flags,:) = trace_temp_array(:,mh_flags,:);

% store result if appropriate
if mcmcInfo.update_flag
    update_index = ceil(mcmcInfo.step/mcmcInfo.update_increment) + 1;
    mcmcInfo.n_steps_inf_array(update_index,:) = mcmcInfo.nStepsCurr;
end