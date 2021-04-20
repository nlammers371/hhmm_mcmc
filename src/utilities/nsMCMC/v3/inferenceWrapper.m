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

    for step = 2:mcmcInfo.n_mcmc_steps %mcmcInfo.n_mcmc_steps          
        waitbar(step/mcmcInfo.n_mcmc_steps,wb);     
        
        mcmcInfo.step = step;

        if ~mcmcInfo.consistencyTestFlag
            % resample chains            
            mcmcInfo = resample_chains_v3(mcmcInfo);    % "Expectation Step"                 

            % perform additional "cross-talk" MH sampling if we are doing chain
            % tempering
            mcmcInfo = temper_chains_v3(mcmcInfo);
        else
            % assign true chains
            sample_chains_true = NaN(size(mcmcInfo.sample_chains));
            for i = 1:length(mcmcInfo.masterSimStruct)
                sample_chains_true(:,:,i) = repmat(mcmcInfo.masterSimStruct(i).naive_states',1,mcmcInfo.n_chains);
            end
            mcmcInfo.sample_chains = sample_chains_true;
        end                
        
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