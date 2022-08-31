function mcmcInfo = initializeVariablesBasicRandom_v2(mcmcInfo)
    
    % initialize transition probabilities
    n_scale = numel(mcmcInfo.observed_fluo)/50;
    n_chains = mcmcInfo.n_chains;
    us_factor = mcmcInfo.upsample_factor;
    
    
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
        
        if ~mcmcInfo.rateSamplingFlag
            
            mcmcInfo.A_curr(:,:,n) = sample_A_dirichlet_par(mcmcInfo.A_alpha(:,:,n), zeros(mcmcInfo.nStates), 1) + 1e-2; % pseuodcounts to avooid underflow issues
            mcmcInfo.A_curr(:,:,n) = mcmcInfo.A_curr(:,:,n) ./ sum(mcmcInfo.A_curr(:,:,n),1);
    %         mcmcInfo.A_curr(:,:,n) = [0.7240    0.2662    0.1379;
    %                                   0.2535    0.6380    0.5591;
    %                                   0.0225    0.0958    0.3030];
        else
            % set hyperparameters
            mcmcInfo.beta_kon = n_scale/2;
            mcmcInfo.beta_koff = n_scale/2;
            mcmcInfo.alpha_kon = 0.1 + rand()*0.25*n_scale/2;
            mcmcInfo.alpha_koff = 0.1 + rand()*0.25*n_scale/2;
            % randomly draw kon and koff
            kon = min([max([0.1 gamrnd(mcmcInfo.alpha_kon,1/mcmcInfo.beta_kon)]) 1]);
            koff = min([max([0.1 gamrnd(mcmcInfo.alpha_koff,1/mcmcInfo.beta_koff)]) 1]);
            % generate rate matrix
            if mcmcInfo.nStates==3
                Q_init = [-2*kon koff 0; 2*kon -kon-koff 2*koff; 0 kon -2*koff];
                mcmcInfo.Q_curr(:,:,n) = Q_init/mcmcInfo.tres*mcmcInfo.upsample_factor;
                % convert to probability matrix (assume single transition)
                mcmcInfo.A_curr(:,:,n) = expm(Q_init);%eye(3) + Q_init + Q_init^2/2 + Q_init^3/6 + Q_init^4/24;
            elseif mcmcInfo.nStates==2
                Q_init = [-kon koff; kon -koff];
                mcmcInfo.Q_curr(:,:,n) = Q_init/mcmcInfo.tres*mcmcInfo.upsample_factor;
                % convert to probability matrix (assume single transition)
                mcmcInfo.A_curr(:,:,n) = eye(2) + Q_init + Q_init^2/2 + Q_init^3/6 + Q_init^4/24;
            else
                error('Rate sampling not supported for nStates>3')
            end
            mcmcInfo.Q_inf_array(:,:,1,n) = mcmcInfo.Q_curr(:,:,n);
        end
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
        mcmcInfo.coeff_MS2_us(:,n) = ms2_loading_coeff_frac(us_factor*mcmcInfo.alphaCurr(n), us_factor*mcmcInfo.nStepsCurr(n), us_factor*mcmcInfo.nStepsMax)';            
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
    