function mcmcInfo = inferenceWrapper(mcmcInfo)

    % initialize chains    
    mcmcInfo = initialize_chains_v4(mcmcInfo);   
    
    % get predicted fluorescence    
    mcmcInfo = predict_fluo_full_v4(mcmcInfo);    
    wb = waitbar(0,'conducting MCMC inference...');

    for step = 2:mcmcInfo.n_mcmc_steps  
        waitbar(step/mcmcInfo.n_mcmc_steps,wb);

        mcmcInfo.step = step;

        % resample chains       
        mcmcInfo = resample_chains_v4(mcmcInfo);    % "Expectation Step"            
        
        % perform additional "cross-talk" MH sampling if we are doing chain
        % tempering
        mcmcInfo = temper_chains_v4(mcmcInfo);
        
        % get updated fluo prediction 
        mcmcInfo = predict_fluo_full_v4(mcmcInfo);
        
        % get empirical transition and occupancy counts    
        mcmcInfo = get_empirical_counts_v4(mcmcInfo);

        % use Gibbs sampling to update hyperparameters      
        mcmcInfo = update_hmm_parameters_v4(mcmcInfo);    

    end
    
    disp('done')
    delete(wb);