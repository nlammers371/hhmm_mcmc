function mcmcInfo = initializeVariablesBasicRandom(mcmcInfo)
    
    % initialize transition probabilities
    n_chains_eff = mcmcInfo.n_chains_eff;
    for n = 1:n_chains_eff    
        alpha_temp = mcmcInfo.A_alpha(:,:,n);
        alpha_temp(eye(mcmcInfo.nStates)==1) = alpha_temp(eye(mcmcInfo.nStates)==1) + normrnd(2,1,mcmcInfo.nStates,1); % distribution hyper params
        mcmcInfo.A_alpha(:,:,n) = 2*alpha_temp;
        mcmcInfo.A_curr(:,:,n) = sample_A_dirichlet_par(mcmcInfo.A_alpha(:,:,n), zeros(mcmcInfo.nStates), 1);
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
    
    for n = 1:n_chains_eff
        mcmcInfo.sigma_curr(n) = sqrt(1./gamrnd(mcmcInfo.a0,1./mcmcInfo.b0));%trandn(-1,Inf)*f_factor/2 + f_factor;%sqrt(1./gamrnd(100*mcmcInfo.seq_length*mcmcInfo.n_traces/2,1./(fluo_vec'*fluo_vec)));
        mcmcInfo.sigma_inf_array(1,n) = mcmcInfo.sigma_curr(n);
    end

    % initialize nSteps
    if ~mcmcInfo.inferNStepsFlag        
      
        mcmcInfo.nStepsCurr = repelem(mcmcInfo.trueParams.nSteps,mcmcInfo.n_chains_eff)'; % current guess (can be fractional)
        mcmcInfo.alphaCurr = mcmcInfo.nStepsCurr * mcmcInfo.alpha_frac;
        
    elseif isfield(mcmcInfo,'nSteps_prior')
      
        mcmcInfo.nStepsCurr = random(mcmcInfo.nSteps_prior,1,mcmcInfo.n_chains)';
        mcmcInfo.alphaCurr = mcmcInfo.nStepsCurr * mcmcInfo.alpha_frac;        
        mcmcInfo.n_steps_inf_array(1,:) = mcmcInfo.nStepsCurr;
        
    else
      
        % generate prior distribution and draw samples        
        mcmcInfo.nStepsCurr = mcmcInfo.nStepsGuess + 2*trandn(-ones(1,mcmcInfo.n_chains_eff),ones(1,mcmcInfo.n_chains_eff));
        mcmcInfo.nStepsCurr(mcmcInfo.nStepsCurr<mcmcInfo.nStepsMin) = mcmcInfo.nStepsMin;
        mcmcInfo.nStepsCurr(mcmcInfo.nStepsCurr>mcmcInfo.nStepsMax) = mcmcInfo.nStepsMax;
        mcmcInfo.alphaCurr = mcmcInfo.nStepsCurr * mcmcInfo.alpha_frac;        
        mcmcInfo.n_steps_inf_array(1,:) = mcmcInfo.nStepsCurr;
        mcmcInfo.nSteps_prior = makedist('Uniform',mcmcInfo.nStepsMin,mcmcInfo.nStepsMax);

    end
    
    % calculate MS2 convolution kernel
    mcmcInfo.coeff_MS2 = NaN(mcmcInfo.nStepsMax,mcmcInfo.n_chains_eff);
    for n = 1:mcmcInfo.n_chains_eff
        mcmcInfo.coeff_MS2(:,n) = ms2_loading_coeff_frac(mcmcInfo.alphaCurr(n), mcmcInfo.nStepsCurr(n), mcmcInfo.nStepsMax)';
    end
    
    % initialize v
    mcmcInfo.v_curr = NaN(mcmcInfo.n_chains_eff,mcmcInfo.nStates);
    v2 = prctile(fluo_vec,95) ./ sum(mcmcInfo.coeff_MS2)';%mean(fluo_vec)/sum(mcmcInfo.coeff_MS2)/(mcmcInfo.pi0_curr(2)+2*mcmcInfo.pi0_curr(3));
    mcmcInfo.v0 = [zeros(n_chains_eff,1) v2];
    if mcmcInfo.nStates==3
         mcmcInfo.v0(:,end+1) = 2*v2;
    end
    mcmcInfo.M0 = eye(mcmcInfo.nStates)*1e1;    
    mcmcInfo.M0(1,1) = 1e5; % add extra weight to "OFF" state
    
    for n = 1:n_chains_eff
        v_cov_mat = 4*mcmcInfo.sigma_curr(n)^2 * inv(mcmcInfo.M0);
        mcmcInfo.v_curr(n,:) = mvnrnd(mcmcInfo.v0(n,:), v_cov_mat)';
        mcmcInfo.v_inf_array(1,:,n) = mcmcInfo.v_curr(n,:);
    end  
    
    