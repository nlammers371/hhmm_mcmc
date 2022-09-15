function mcmcInfo = ga_move(mcmcInfo)    
    step = mcmcInfo.step;
    logL_vec = mcmcInfo.logL_vec(step-1,:)*numel(mcmcInfo.observed_fluo);

%             rs_factor = exp(2*logsumexp(logL_vec,2) - logsumexp(2*logL_vec,2));
    if true%rs_factor <= mcmcInfo.n_chains/2
        logL_vec = logL_vec-logsumexp(logL_vec,2);
%             [~,si_vec] = sort(logL_vec,'descend');
        prob_vec = exp(logL_vec); % inject pseudocounts to avoid sparsity issues
%             rs_ids = repelem(si_vec(1:5),5);              
        rs_ids = randsample(1:mcmcInfo.n_chains,mcmcInfo.n_chains,true,prob_vec);          
        % resample latent states and fluo predictions
        mcmcInfo.sample_chains = mcmcInfo.sample_chains(:,rs_ids,:);
        mcmcInfo.sample_fluo = mcmcInfo.sample_fluo(:,rs_ids,:);
        mcmcInfo.A_curr = mcmcInfo.A_curr(:,:,rs_ids);
        if mcmcInfo.rateSamplingFlag
            mcmcInfo.Q_curr = mcmcInfo.Q_curr(:,:,rs_ids);
        end
        mcmcInfo.v_curr = mcmcInfo.v_curr(rs_ids,:);
        mcmcInfo.pi0_curr = mcmcInfo.pi0_curr(rs_ids,:);
        mcmcInfo.sigma_curr = mcmcInfo.sigma_curr(rs_ids);
    end