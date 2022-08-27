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
        
        if mcmcInfo.resampleTracesFlag && step >= mcmcInfo.burn_in && mod(step,mcmcInfo.rs_freq) == 0
            mcmcInfo = ga_move(mcmcInfo);
        end
        
        % if tempering, artificially raise sigma values
        if mcmcInfo.annealingSigmaFlag
            mcmcInfo.sigma_curr = mcmcInfo.sigma_curr .* mcmcInfo.tempGradArray(mcmcInfo.step,:)';
        end
        
        % resample chains  promoter state sequences        
        if mcmcInfo.mhResamplingFlag
            mcmcInfo = resample_chains_v5_MH(mcmcInfo);
        elseif false%mcmcInfo.rateSamplingFlag
            mcmcInfo = resample_chains_v5_MH(mcmcInfo);                 
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