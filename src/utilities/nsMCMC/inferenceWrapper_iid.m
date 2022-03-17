function mcmcInfo = inferenceWrapper_iid(mcmcInfo)
    
    % initialize step
    mcmcInfo.step = 1;
    
    % initialize chains    
    mcmcInfo = initialize_chains_iid(mcmcInfo);    

    % get predicted fluorescence
    mcmcInfo = predict_fluo_full_v3(mcmcInfo);

    % calculate initial logL
    mcmcInfo = calculateLogLikelihood_iid(mcmcInfo);
    
    wb = waitbar(0,'conducting iid MCMC inference...');

    for step = 2:mcmcInfo.n_mcmc_steps    
        waitbar(step/mcmcInfo.n_mcmc_steps,wb);     
        
        mcmcInfo.step = step;       
          
        % resample chains  promoter state sequences        
        mcmcInfo = resample_chains_iid(mcmcInfo);    % "Expectation Step"                 
                       
        % get empirical transition and occupancy counts    
        mcmcInfo = get_empirical_counts_iid(mcmcInfo);

        % use Gibbs sampling to update hyperparameters      
        mcmcInfo = update_hmm_parameters_iid(mcmcInfo);                
        
        % calculate updated logL
        mcmcInfo = calculateLogLikelihood_iid(mcmcInfo);

    end
    
    disp('done')
    delete(wb);