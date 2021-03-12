function mcmcInfo = inferenceWrapper(mcmcInfo)

    % initialize chains
    mcmcInfo = initialize_chains_v3(mcmcInfo);

    % get predicted fluorescence
    mcmcInfo = predict_fluo_full(mcmcInfo);

    wb = waitbar(0,'conducting MCMC inference...');

    for step = 2:mcmcInfo.n_mcmc_steps %mcmcInfo.n_mcmc_steps    
        waitbar(step/mcmcInfo.n_mcmc_steps,wb);

        mcmcInfo.step = step;

        % resample chains    
        if mcmcInfo.MHResampling
            mcmcInfo = resample_chains_MH(mcmcInfo);              
        else
            mcmcInfo = resample_chains_v3(mcmcInfo);    % "Expectation Step"         
        end
        
        % get empirical transition and occupancy counts    
        mcmcInfo = get_empirical_counts_v3(mcmcInfo);

        % use Gibbs sampling to update hyperparameters      
        mcmcInfo = update_hmm_parameters_v3(mcmcInfo);    

    end
    
    disp('done')
    delete(wb);