function mcmcInfo = intitializeVariablesTrue(mcmcInfo)
    
    % A prior--assume strongly diagonal PDF given short timescale
    % take A columns to follow multinomial Dirichlet distribution
    mcmcInfo.A_curr = NaN(mcmcInfo.nStates,mcmcInfo.nStates,mcmcInfo.n_chains);
    mcmcInfo.pi0_curr = NaN(mcmcInfo.n_chains,mcmcInfo.nStates);
    mcmcInfo.A_alpha = ones(mcmcInfo.nStates,mcmcInfo.nStates,mcmcInfo.n_chains);%*n_particles*n_traces;

    n_chains = mcmcInfo.n_chains;
    for n = 1:n_chains    
        alpha_temp = mcmcInfo.A_alpha(:,:,n);
        alpha_temp(eye(mcmcInfo.nStates)==1) = alpha_temp(eye(mcmcInfo.nStates)==1) + rand(mcmcInfo.nStates,1)*10; % distribution hyper params
        mcmcInfo.A_alpha(:,:,n) = alpha_temp;
        mcmcInfo.A_curr(:,:,n) = mcmcInfo.trueParams.A;
        mcmcInfo.A_inf_array(:,:,1,n) = mcmcInfo.A_curr(:,:,n);

        % calculate pi0 
        [V, D] = eig(mcmcInfo.A_curr(:,:,n));
        [~, mi] = max(real(diag(D)));
        mcmcInfo.pi0_curr(n,:) = V(:,mi)/sum(V(:,mi));
        mcmcInfo.pi0_inf_array(1,:,n) = mcmcInfo.pi0_curr(n,:);
    end

    % initialize sigma as inverse gamma (see: http://ljwolf.org/teaching/gibbs.html)
    for n = 1:n_chains
        mcmcInfo.sigma_curr(n) = mcmcInfo.trueParams.sigma;%sqrt(1./gamrnd(100*mcmcInfo.seq_length*mcmcInfo.n_traces/2,1./(fluo_vec'*fluo_vec)));
        mcmcInfo.sigma_inf_array(1,n) = mcmcInfo.sigma_curr(n);
    end

    % initialize v
    mcmcInfo.v_curr = NaN(mcmcInfo.n_chains,mcmcInfo.nStates);    
    for n = 1:n_chains
        mcmcInfo.v_curr(n,:) = mcmcInfo.trueParams.v;
        mcmcInfo.v_inf_array(1,:,n) = mcmcInfo.v_curr(n,:);
    end