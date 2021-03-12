function mcmcInfo = update_hmm_parameters(mcmcInfo)    

    % extrace parameters
    nStates = mcmcInfo.nStates;
    seq_length = mcmcInfo.seq_length;
    n_traces = mcmcInfo.n_traces;
%     n_chains = mcmcInfo.n_chains;
    coeff_MS2 = mcmcInfo.coeff_MS2;
    
    % update A
    A_counts = sum(mcmcInfo.transition_count_array,3)/mcmcInfo.n_chains;    
    mcmcInfo.A_curr = sample_A_dirichlet(mcmcInfo.A_alpha, A_counts);
    mcmcInfo.A_inf_array(:,:,mcmcInfo.step) = mcmcInfo.A_curr;
    
    % update V               
         
    % get current likelihood
    v_orig = mcmcInfo.v_curr;
    logL_fluo_orig = calculate_fluo_logL_full(mcmcInfo);

    % draw sample
    v_new = mcmcInfo.v_curr;
    prop_sigma = max([0.1*v_new,repmat(0.025,nStates,1)],[],2);
    v_new = normrnd(v_new,prop_sigma);%trandn(-v_new(n),Inf)*mcmcInfo.v_prop_sigma + v_new(n);

    % calculate updated likelihood
    mcmcInfo.v_curr = v_new;
    logL_fluo_new = calculate_fluo_logL_full(mcmcInfo);

    % perform MH move
    accept_flag = exp(logL_fluo_new-logL_fluo_orig) > rand();

    % reset value if move not accepted
    if ~accept_flag
        mcmcInfo.v_curr = v_orig;
    end

    % record outcome
    mcmcInfo.v_acceptance_array(mcmcInfo.step) = accept_flag;
   
    mcmcInfo.v_inf_array(mcmcInfo.step,:) = mcmcInfo.v_curr;
    
    % fix sigma for now
    mcmcInfo.sigma_curr = mcmcInfo.sigma;
    