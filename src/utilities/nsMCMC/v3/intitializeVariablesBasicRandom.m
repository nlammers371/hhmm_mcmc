function mcmcInfo = intitializeVariablesBasicRandom(mcmcInfo)
    
    % A prior--assume strongly diagonal PDF given short timescale
    % take A columns to follow multinomial Dirichlet distribution
    mcmcInfo.A_curr = NaN(mcmcInfo.nStates,mcmcInfo.nStates,mcmcInfo.n_chains*mcmcInfo.n_temps_per_chain);
    mcmcInfo.pi0_curr = NaN(mcmcInfo.n_chains*mcmcInfo.n_temps_per_chain,mcmcInfo.nStates);
    mcmcInfo.A_alpha = ones(mcmcInfo.nStates,mcmcInfo.nStates,mcmcInfo.n_chains*mcmcInfo.n_temps_per_chain);%*n_particles*n_traces;

    n_chains = mcmcInfo.n_chains_eff;
    for n = 1:n_chains    
        alpha_temp = mcmcInfo.A_alpha(:,:,n);
        alpha_temp(eye(mcmcInfo.nStates)==1) = alpha_temp(eye(mcmcInfo.nStates)==1) + rand(mcmcInfo.nStates,1)*10; % distribution hyper params
        mcmcInfo.A_alpha(:,:,n) = alpha_temp;
        mcmcInfo.A_curr(:,:,n) = sample_A_dirichlet_par(alpha_temp, zeros(mcmcInfo.nStates), 1);
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
        mcmcInfo.sigma_curr(n) = sqrt(1./gamrnd(mcmcInfo.a0,1./mcmcInfo.b0));%trandn(-1,Inf)*f_factor/2 + f_factor;%sqrt(1./gamrnd(100*mcmcInfo.seq_length*mcmcInfo.n_traces/2,1./(fluo_vec'*fluo_vec)));
        mcmcInfo.sigma_inf_array(1,n) = mcmcInfo.sigma_curr(n);
    end

    % initialize v
    mcmcInfo.v_curr = NaN(mcmcInfo.n_chains*mcmcInfo.n_temps_per_chain,mcmcInfo.nStates);
    v2 = prctile(fluo_vec,99) / mcmcInfo.nSteps;%mean(fluo_vec)/sum(mcmcInfo.coeff_MS2)/(mcmcInfo.pi0_curr(2)+2*mcmcInfo.pi0_curr(3));
    mcmcInfo.v0 = [0 v2 2*v2];
    mcmcInfo.M0 = eye(mcmcInfo.nStates)*1e1;    
    mcmcInfo.M0(1,1) = 1e5;
    
    for n = 1:n_chains
        v_cov_mat = mcmcInfo.sigma_curr(n)^2 * inv(mcmcInfo.M0);
        mcmcInfo.v_curr(n,:) = mvnrnd(mcmcInfo.v0, v_cov_mat)';
        mcmcInfo.v_inf_array(1,:,n) = mcmcInfo.v_curr(n,:);
    end  