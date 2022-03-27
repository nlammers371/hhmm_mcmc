function mcmcInfo = initializeVariablesBasicRandom(mcmcInfo)
    
    % initialize transition probabilities
%     mcmcInfo.A_alpha = ones(mcmcInfo.nStates,mcmcInfo.nStates,mcmcInfo.n_chains_eff);
    n_scale = numel(mcmcInfo.observed_fluo)/50;
    n_chains_eff = mcmcInfo.n_chains_eff;
    for n = 1:n_chains_eff    
%         alpha_temp = mcmcInfo.A_alpha(:,:,n);        
        alpha_temp = abs(normrnd(n_scale,n_scale/3,mcmcInfo.nStates,mcmcInfo.nStates));        
        alpha_temp(eye(mcmcInfo.nStates)~=1) = alpha_temp(eye(mcmcInfo.nStates)~=1)/2; % distribution hyper params
        if mcmcInfo.nStates == 3
            alpha_temp(1,3) = alpha_temp(1,3)/4;
            alpha_temp(3,1) = alpha_temp(3,1)/4;
        end
        if ~mcmcInfo.strongAPriorFlag
            alpha_temp = alpha_temp / 2;
        end
        mcmcInfo.A_alpha(:,:,n) = alpha_temp;
        mcmcInfo.A_curr(:,:,n) = sample_A_dirichlet_par(mcmcInfo.A_alpha(:,:,n), zeros(mcmcInfo.nStates), 1);
        mcmcInfo.A_inf_array(:,:,1,n) = mcmcInfo.A_curr(:,:,n);
        
        % calculate pi0 
        [V, D] = eig(mcmcInfo.A_curr(:,:,n));
        [~, mi] = max(real(diag(D)));
        mcmcInfo.pi0_curr(n,:) = V(:,mi)/sum(V(:,mi));
        mcmcInfo.pi0_inf_array(1,:,n) = mcmcInfo.pi0_curr(n,:);
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
    fluo_vec_mid = fluo_vec(fluo_vec >= mean(fluo_vec) & fluo_vec <= mean(fluo_vec) + 2*std(fluo_vec));
    v2 = randsample(fluo_vec_mid,n_chains_eff,true) ./ sum(mcmcInfo.coeff_MS2)';%mean(fluo_vec)/sum(mcmcInfo.coeff_MS2)/(mcmcInfo.pi0_curr(2)+2*mcmcInfo.pi0_curr(3));
    mcmcInfo.v0 = [zeros(n_chains_eff,1) v2];
    if mcmcInfo.nStates==3
         mcmcInfo.v0(:,end+1) = 2*v2;
    end
    mcmcInfo.M0 = eye(mcmcInfo.nStates)*1;    
    mcmcInfo.M0(1,1) = 1e5; % add extra weight to "OFF" state
    
    for n = 1:n_chains_eff
        v_cov_mat = mcmcInfo.sigma_curr(n)^2 * inv(mcmcInfo.M0);
        mcmcInfo.v_curr(n,:) = sort(mvnrnd(mcmcInfo.v0(n,:), v_cov_mat)');
        mcmcInfo.v_inf_array(1,:,n) = mcmcInfo.v_curr(n,:);
    end      
    
    % if we're doing ensemble inference, reset everything to have single
    % value
%     if mcmcInfo.ensembleInferenceFlag
%         mcmcInfo.v0 = repmat(mcmcInfo.v0(1,:),mcmcInfo.n_chains_eff,1);
%         mcmcInfo.v_curr = repmat(mcmcInfo.v_curr(1,:),mcmcInfo.n_chains_eff,1);
%         mcmcInfo.pi0_curr = repmat(mcmcInfo.pi0_curr(1,:),mcmcInfo.n_chains_eff,1);
%         mcmcInfo.sigma_curr = repmat(mcmcInfo.sigma_curr(1),mcmcInfo.n_chains_eff,1);
%         mcmcInfo.A_curr = repmat(mcmcInfo.A_curr(:,:,1),1,1,mcmcInfo.n_chains_eff);
%         mcmcInfo.A_alpha = repmat(mcmcInfo.A_alpha(:,:,1),1,1,mcmcInfo.n_chains_eff);
%         mcmcInfo.coeff_MS2 = repmat(mcmcInfo.coeff_MS2(:,1),1,mcmcInfo.n_chains_eff);
%         mcmcInfo.nStepsCurr = repmat(mcmcInfo.nStepsCurr(1),mcmcInfo.n_chains_eff,1);
%         mcmcInfo.alphaCurr = repmat(mcmcInfo.alphaCurr(1),mcmcInfo.n_chains_eff,1);
%     end   
        