function mcmcInfo = inferenceWrapper_PF(mcmcInfo)
    
    % initialize step
    mcmcInfo.step = 1;
    
    % initialize chains
    if ~isfield(mcmcInfo, 'sample_chains')
        mcmcInfo = initialize_chains_v3(mcmcInfo);
    end

    % get predicted fluorescence
    mcmcInfo = predict_fluo_full_v3(mcmcInfo);

    % calculate initial logL
    mcmcInfo = calculateLogLikelihood(mcmcInfo);
    
    wb = waitbar(0,'conducting MCMC inference...');

    for step = 2:mcmcInfo.n_mcmc_steps    
        waitbar(step/mcmcInfo.n_mcmc_steps,wb);     
        
        mcmcInfo.step = step;       
        
        % if tempering, artificially raise sigma values
        if mcmcInfo.annealingSigmaFlag
            mcmcInfo.sigma_curr = mcmcInfo.sigma_curr .* mcmcInfo.tempGradArray(mcmcInfo.step,:)';
        end
        
        % resample chains  promoter state sequences        
        mcmcInfo = resample_chains_PF(mcmcInfo);

        % get empirical transition and occupancy counts    
        mcmcInfo = get_empirical_counts_v3(mcmcInfo);

        % use Gibbs sampling to update hyperparameters           
        mcmcInfo = update_hmm_parameters_v3(mcmcInfo);   
  
        % if desired, perform MH moves to sample "nSteps"
        if mcmcInfo.inferNStepsFlag
            mcmcInfo = mh_sample_nSteps(mcmcInfo);           
        end
        
        % calculate updated logL
        mcmcInfo = calculateLogLikelihood(mcmcInfo);

    end
    
    disp('done')
    delete(wb);