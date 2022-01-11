function mcmcInfo = initializeVariablesWithPriors(mcmcInfo,mcmcInfoInit,ids_to_use)
    
    % initialize transition probabilities
    n_chains = mcmcInfo.n_chains_eff;
    for n = 1:n_chains    
        mcmcInfo.A_alpha(:,:,n) = mcmcInfoInit.A_alpha(:,:,ids_to_use(n));
        mcmcInfo.A_curr(:,:,n) = mcmcInfoInit.A_curr(:,:,ids_to_use(n));
        mcmcInfo.A_inf_array(:,:,1,n) = mcmcInfo.A_curr(:,:,n);

        % calculate pi0 
        [V, D] = eig(mcmcInfo.A_curr(:,:,n));
        [~, mi] = max(real(diag(D)));
        mcmcInfo.pi0_curr(n,:) = V(:,mi)/sum(V(:,mi));
        mcmcInfo.pi0_inf_array(1,:,n) = mcmcInfo.pi0_curr(n,:);
    end

    % initialize sigma as inverse gamma (see: http://ljwolf.org/teaching/gibbs.html)
    mcmcInfo.sigma_curr = NaN(mcmcInfo.n_chains*mcmcInfo.n_temps_per_chain,1);
    mcmcInfo.a0 = numel(mcmcInfo.observed_fluo)/100;
    fluo_vec = mcmcInfo.observed_fluo(:);
    mcmcInfo.b0 = 0.3*mean(fluo_vec).^2;
    
    for n = 1:n_chains
        mcmcInfo.sigma_curr(n) = mcmcInfoInit.sigma_curr(ids_to_use(n));
        mcmcInfo.sigma_inf_array(1,n) = mcmcInfo.sigma_curr(n);
    end

    % initialize nSteps
    if ~mcmcInfo.inferNStepsFlag
        mcmcInfo.nStepsCurr = repelem(mcmcInfo.trueParams.nSteps,mcmcInfo.n_chains_eff)'; % current guess (can be fractional)
        mcmcInfo.alphaCurr = mcmcInfo.nStepsCurr * mcmcInfo.alpha_frac;
    else
        % generate prior distribution and draw samples        
        mcmcInfo.nStepsCurr = mcmcInfoInit.nStepsCurr(ids_to_use);
        mcmcInfo.alphaCurr = mcmcInfo.nStepsCurr * mcmcInfo.alpha_frac;
        mcmcInfo.n_steps_inf_array(1,:) = mcmcInfo.nStepsCurr;
    end
    
    % calculate MS2 convolution kernel
    mcmcInfo.coeff_MS2  = NaN(mcmcInfo.nStepsMax,mcmcInfo.n_chains_eff);
    for n = 1:mcmcInfo.n_chains_eff
        mcmcInfo.coeff_MS2(:,n) = ms2_loading_coeff_frac(mcmcInfo.alphaCurr(n), mcmcInfo.nStepsCurr(n), mcmcInfo.nStepsMax)';
    end
    
    % initialize v
    mcmcInfo.v_curr = NaN(mcmcInfo.n_chains_eff,mcmcInfo.nStates);    
    mcmcInfo.v0 = mcmcInfo.v0(ids_to_use,:);
    mcmcInfo.M0 = eye(mcmcInfo.nStates)*1e1;    
    mcmcInfo.M0(1,1) = 1e5; % add extra weight to "OFF" state
    
    for n = 1:n_chains        
        mcmcInfo.v_curr(n,:) = mcmcInfoInit.v_curr(ids_to_use(n),:);
        mcmcInfo.v_inf_array(1,:,n) = mcmcInfo.v_curr(n,:);
    end  
    
    