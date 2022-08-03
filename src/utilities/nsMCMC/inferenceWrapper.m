function mcmcInfo = inferenceWrapper(mcmcInfo)
    
    % initialize step
    mcmcInfo.step = 1;
    
    % initialize chains    
    mcmcInfo = initialize_chains_v3(mcmcInfo);    

    % get predicted fluorescence
    mcmcInfo = predict_fluo_full_v3(mcmcInfo);

    % calculate initial logL
    mcmcInfo = calculateLogLikelihood(mcmcInfo);
    
    wb = waitbar(0,'conducting MCMC inference...');

    for step = 2:mcmcInfo.n_mcmc_steps    
        waitbar(step/mcmcInfo.n_mcmc_steps,wb);     
        mcmcInfo.step = step;       
        
        if mcmcInfo.resampleTracesFlag && step >= mcmcInfo.burn_in %&& mod(step,mcmcInfo.rs_freq) == 0
            % resample using total likelihood
            logL_vec = mcmcInfo.logL_vec(step-1,:)*numel(mcmcInfo.observed_fluo);
            
            rs_factor = exp(2*logsumexp(logL_vec,2) - logsumexp(2*logL_vec,2));
            if true%rs_factor <= mcmcInfo.n_chains/2
                logL_vec = logL_vec-logsumexp(logL_vec,2);
%             [~,si_vec] = sort(logL_vec,'descend');
                prob_vec = exp(logL_vec); % inject pseudocounts to avoid sparsity issues
    %             rs_ids = repelem(si_vec(1:5),5);
                try
                    rs_ids = randsample(1:mcmcInfo.n_chains,mcmcInfo.n_chains,true,prob_vec);
                catch
                    error('check')
                end
                % resample latent states and fluo predictions
                mcmcInfo.sample_chains = mcmcInfo.sample_chains(:,rs_ids,:);
                mcmcInfo.sample_fluo = mcmcInfo.sample_fluo(:,rs_ids,:);
                mcmcInfo.A_curr = mcmcInfo.A_curr(:,:,rs_ids);
                mcmcInfo.v_curr = mcmcInfo.v_curr(rs_ids,:);
                mcmcInfo.pi0_curr = mcmcInfo.pi0_curr(rs_ids,:);
                mcmcInfo.sigma_curr = mcmcInfo.sigma_curr(rs_ids);
            end
        end
        
        % if tempering, artificially raise sigma values
        if mcmcInfo.annealingSigmaFlag
            mcmcInfo.sigma_curr = mcmcInfo.sigma_curr .* mcmcInfo.tempGradArray(mcmcInfo.step,:)';
        end
        
        % resample chains  promoter state sequences        
        if mcmcInfo.mhResamplingFlag
            mcmcInfo = resample_chains_v5_MH(mcmcInfo);
        elseif mcmcInfo.PFResamplingFlag
            mcmcInfo = resample_chains_PF_v3(mcmcInfo);                 
        else
            mcmcInfo = resample_chains_v5(mcmcInfo);
        end
        
        % get empirical transition and occupancy counts          
        mcmcInfo = get_empirical_counts_v3(mcmcInfo);   
        
        % use Gibbs sampling to update hyperparameters           
        mcmcInfo = update_hmm_parameters_v4(mcmcInfo);          
        
        % if desired, perform MH moves to sample "nSteps"
        if mcmcInfo.inferNStepsFlag
            mcmcInfo = mh_sample_nSteps(mcmcInfo);           
        end
        
        % calculate updated logL        
        mcmcInfo = calculateLogLikelihood(mcmcInfo);        
%         toc
    end
    
    disp('done')
    delete(wb);