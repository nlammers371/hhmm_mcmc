function mcmcInfo = initialize_mcmc_parameters(mcmcInfo)

    % initialize sigma as inverse gamma (see: http://ljwolf.org/teaching/gibbs.html)
    fluo_vec = mcmcInfo.observed_fluo(:);
    f_sigma = 0.05*mean(fluo_vec);
    if ~mcmcInfo.par_chain_flag
        mcmcInfo.sigma_curr = trandn(-2,Inf)*f_sigma + 2*f_sigma;%sqrt(1./gamrnd(100*mcmcInfo.seq_length*mcmcInfo.n_traces/2,1./(fluo_vec'*fluo_vec)));
    else
        mcmcInfo.sigma_curr = trandn(repmat(-2,mcmcInfo.n_chains,1),Inf(mcmcInfo.n_chains,1))*f_sigma + 2*f_sigma;
    end

    %% initialize v
    v2 = prctile(fluo_vec,96) / mcmcInfo.nSteps;%mean(fluo_vec)/sum(mcmcInfo.coeff_MS2)/(mcmcInfo.pi0_curr(2)+2*mcmcInfo.pi0_curr(3));
    v2_err = std(fluo_vec) / mcmcInfo.nSteps;
    v_vec = [0 ; v2 ; 2*v2];
    if ~mcmcInfo.par_chain_flag
        mcmcInfo.v_curr = [0 ; v2 ; 2*v2] + rand(3,1)*.2;
    else
        mcmcInfo.v_curr = zeros(mcmcInfo.n_chains,mcmcInfo.nStates);
        for k = 2:mcmcInfo.nStates
            mcmcInfo.v_curr(:,k) = trandn(repmat(-v_vec(k)/v2_err,mcmcInfo.n_chains,1) + mcmcInfo.v_curr(:,k-1)/v2_err,Inf(mcmcInfo.n_chains,1))*v2_err + v_vec(k);
        end
    end
    % A prior--assume strongly diagonal PDF given short timescale
    % take A columns to follow multinomial Dirichlet distribution
    mcmcInfo.A_alpha = ones(mcmcInfo.nStates);%*n_particles*n_traces;
    mcmcInfo.A_alpha(eye(mcmcInfo.nStates)==1) = mcmcInfo.A_alpha(eye(mcmcInfo.nStates)==1)*10; % distribution hyper params
    if ~mcmcInfo.par_chain_flag
        mcmcInfo.A_curr = sample_A_dirichlet(mcmcInfo.A_alpha, zeros(mcmcInfo.nStates),1);
    else
        mcmcInfo.A_curr = sample_A_dirichlet(mcmcInfo.A_alpha, zeros(mcmcInfo.nStates),mcmcInfo.n_chains);
    end
    % calculate pi0
    if ~mcmcInfo.par_chain_flag 
        [V, D] = eig(mcmcInfo.A_curr);
        [~, mi] = max(real(diag(D)));
        mcmcInfo.pi0_curr = V(:,mi)/sum(V(:,mi));
    else
        mcmcInfo.pi0_curr = NaN(mcmcInfo.n_chains,mcmcInfo.nStates);
        for n = 1:mcmcInfo.n_chains
            [V, D] = eig(mcmcInfo.A_curr(:,:,n));
            [~, mi] = max(real(diag(D)));
            mcmcInfo.pi0_curr(n,:) = V(:,mi)/sum(V(:,mi));
        end
    end