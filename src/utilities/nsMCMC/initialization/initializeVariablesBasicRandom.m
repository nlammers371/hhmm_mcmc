function mcmcInfo = initializeVariablesBasicRandom(mcmcInfo)
    
    % initialize transition probabilities
    n_scale = numel(mcmcInfo.observed_fluo)/50;
    n_chains = mcmcInfo.n_chains;
    
    if ~mcmcInfo.reducedModelFlag && ~mcmcInfo.mhInferenceFlag
        for n = 1:n_chains    
            alpha_temp = abs(normrnd(n_scale,n_scale/3,mcmcInfo.nStates,mcmcInfo.nStates));        

            if mcmcInfo.nStates == 3
                alpha_temp(1,3) = alpha_temp(1,3)/4;
                alpha_temp(3,1) = alpha_temp(3,1)/4;
            end
            if ~mcmcInfo.strongAPriorFlag
                alpha_temp = alpha_temp / 2;
            end
            mcmcInfo.A_alpha(:,:,n) = alpha_temp;
            mcmcInfo.A_curr(:,:,n) = sample_A_dirichlet_par(mcmcInfo.A_alpha(:,:,n), zeros(mcmcInfo.nStates), 1) + 1e-2; % pseuodcounts to avooid underflow issues
            mcmcInfo.A_curr(:,:,n) = mcmcInfo.A_curr(:,:,n) ./ sum(mcmcInfo.A_curr(:,:,n),1);
    %         mcmcInfo.A_curr(:,:,n) = [0.7240    0.2662    0.1379;
    %                                   0.2535    0.6380    0.5591;
    %                                   0.0225    0.0958    0.3030];
            mcmcInfo.A_inf_array(:,:,1,n) = mcmcInfo.A_curr(:,:,n);

            % calculate pi0 
            [V, D] = eig(mcmcInfo.A_curr(:,:,n));
            [~, mi] = max(real(diag(D)));
            mcmcInfo.pi0_curr(n,:) = V(:,mi)/sum(V(:,mi));
            mcmcInfo.pi0_inf_array(1,:,n) = mcmcInfo.pi0_curr(n,:);

        end

        % initialize sigma as inverse gamma (see: http://ljwolf.org/teaching/gibbs.html)    
        mcmcInfo.a0 = numel(mcmcInfo.observed_fluo)/2;
        fluo_vec = mcmcInfo.observed_fluo(:);
        mcmcInfo.b0 = numel(mcmcInfo.observed_fluo)*(0.25*mean(fluo_vec)).^2;

        for n = 1:n_chains
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
        v2 = randsample(fluo_vec_mid,n_chains,true) ./ sum(mcmcInfo.coeff_MS2)';
        mcmcInfo.v0 = [zeros(n_chains,1) v2];
        if mcmcInfo.nStates==3
             mcmcInfo.v0(:,end+1) = 2.5*v2;
        elseif mcmcInfo.nStates==2
            mcmcInfo.v0(:,2) = 2*v2;
        end
        mcmcInfo.M0 = eye(mcmcInfo.nStates)*1;    
        mcmcInfo.M0(1,1) = 1e5; % add extra weight to "OFF" state

        for n = 1:n_chains
            v_cov_mat = mcmcInfo.sigma_curr(n)^2 * inv(mcmcInfo.M0);
            v_temp = sort(mvnrnd(mcmcInfo.v0(n,:), v_cov_mat)');        
            v_temp(v_temp<0) = 0;

            mcmcInfo.v_curr(n,:) = v_temp; %[0 2 5];%
            mcmcInfo.v_inf_array(1,:,n) = mcmcInfo.v_curr(n,:);
        end      
    elseif mcmcInfo.reducedModelFlag && mcmcInfo.mhInferenceFlag && mcmcInfo.nStates == 3
        kon_init = rand()*0.02;     
        koff_init = rand()*0.02;     
        k_corr_init = rand()-0.5;
            
        mcmcInfo.k_curr = [kon_init,koff_init,k_corr_init];
        
        Q_init = Q_helper_fun(kon_init,koff_init,exp(k_corr_init));
               
        mcmcInfo.A_curr = repmat(expm(Q_init*mcmcInfo.tres),1,1,n_chains);

        mcmcInfo.k_inf_array(1,:) = [kon_init koff_init k_corr_init];

        % calculate pi0 
        [V, D] = eig(mcmcInfo.A_curr(:,:,1));
        [~, mi] = max(real(diag(D)));
        mcmcInfo.pi0_curr = repmat((V(:,mi)/sum(V(:,mi)))',n_chains,1);                                        
        
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
        fluo_vec = mcmcInfo.observed_fluo(:);
        fluo_vec_mid = fluo_vec(fluo_vec >= mean(fluo_vec) & fluo_vec <= mean(fluo_vec) + 2*std(fluo_vec));
        v2 = randsample(fluo_vec_mid,1,true) ./ sum(mcmcInfo.coeff_MS2(:,1))';
        v_corr_init = rand()-0.5;
        vc = exp(v_corr_init);
        mcmcInfo.v_curr = repmat([0 v2 2*v2*vc],n_chains,1);
        
        % fluorescence noise parameter
        mcmcInfo.sigma_curr = repmat(0.2*mean(fluo_vec),n_chains,1);
        
        mcmcInfo.r_curr = [v2/mcmcInfo.tres v_corr_init];
        mcmcInfo.r_inf_array(1,:) = mcmcInfo.r_curr;
    else
        error('Incompatible inference options')
    end