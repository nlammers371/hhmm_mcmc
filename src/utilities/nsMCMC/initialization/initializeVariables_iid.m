function mcmcInfo = initializeVariables_iid(mcmcInfo)
    
    % initialize transition probabilities
    mcmcInfo.pi0_alpha = NaN(mcmcInfo.n_chains_eff,mcmcInfo.nStates);
    n_chains_eff = mcmcInfo.n_chains_eff;
    for n = 1:n_chains_eff    
        pi0_temp = ones(1,mcmcInfo.nStates);
        pi0_samp = drchrnd(pi0_temp,1);             
        mcmcInfo.pi0_alpha(n,:) = pi0_samp*0.5*mcmcInfo.n_traces*mcmcInfo.seq_length;   
        mcmcInfo.pi0_curr(n,:) = pi0_samp;
    end

    % initialize sigma as inverse gamma (see: http://ljwolf.org/teaching/gibbs.html)    
    mcmcInfo.a0 = numel(mcmcInfo.observed_fluo)/100;
    fluo_vec = mcmcInfo.observed_fluo(:);
    mcmcInfo.b0 = 0.3*mean(fluo_vec).^2;
    
    for n = 1:n_chains_eff
        mcmcInfo.sigma_curr(n) = sqrt(1./gamrnd(mcmcInfo.a0,1./mcmcInfo.b0));%trandn(-1,Inf)*f_factor/2 + f_factor;%sqrt(1./gamrnd(100*mcmcInfo.seq_length*mcmcInfo.n_traces/2,1./(fluo_vec'*fluo_vec)));
        mcmcInfo.sigma_inf_array(1,n) = mcmcInfo.sigma_curr(n);
    end

    % initialize nSteps
    if ~mcmcInfo.inferNStepsFlag
        mcmcInfo.nStepsCurr = repelem(mcmcInfo.nSteps,mcmcInfo.n_chains_eff)'; % current guess (can be fractional)
        mcmcInfo.alphaCurr = mcmcInfo.nStepsCurr * mcmcInfo.alpha_frac;
    else
        % generate prior distribution and draw samples        
        mcmcInfo.nStepsCurr = mcmcInfo.nStepsGuess + 2*trandn(-ones(1,mcmcInfo.n_chains_eff),ones(1,mcmcInfo.n_chains_eff));
        mcmcInfo.nStepsCurr(mcmcInfo.nStepsCurr<mcmcInfo.nStepsMin) = mcmcInfo.nStepsMin;
        mcmcInfo.nStepsCurr(mcmcInfo.nStepsCurr>mcmcInfo.nStepsMax) = mcmcInfo.nStepsMax;
        mcmcInfo.alphaCurr = mcmcInfo.nStepsCurr * mcmcInfo.alpha_frac;
        mcmcInfo.n_steps_inf_array(1,:) = mcmcInfo.nStepsCurr;
    end
    
    % calculate MS2 convolution kernel
    mcmcInfo.coeff_MS2 = NaN(mcmcInfo.nStepsMax,mcmcInfo.n_chains_eff);
    for n = 1:mcmcInfo.n_chains_eff
        mcmcInfo.coeff_MS2(:,n) = ms2_loading_coeff_frac(mcmcInfo.alphaCurr(n), mcmcInfo.nStepsCurr(n), mcmcInfo.nStepsMax)';
    end
    
    % initialize v    
    v2 = prctile(fluo_vec,99) ./ sum(mcmcInfo.coeff_MS2)';%mean(fluo_vec)/sum(mcmcInfo.coeff_MS2)/(mcmcInfo.pi0_curr(2)+2*mcmcInfo.pi0_curr(3));
    mcmcInfo.v0 = [zeros(n_chains_eff,1) v2];
    if mcmcInfo.nStates==3
         mcmcInfo.v0(:,end+1) = 2*v2;
    end
    mcmcInfo.M0 = eye(mcmcInfo.nStates)*0.4;    
    mcmcInfo.M0(1,1) = 1e5; % add extra weight to "OFF" state
    
    for n = 1:n_chains_eff
        v_cov_mat = mcmcInfo.sigma_curr(n)^2 * inv(mcmcInfo.M0);
        mcmcInfo.v_curr(n,:) = mvnrnd(mcmcInfo.v0(n,:), v_cov_mat)';
        mcmcInfo.v_inf_array(1,:,n) = mcmcInfo.v_curr(n,:);
    end         