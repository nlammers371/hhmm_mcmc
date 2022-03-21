function mcmcInfo = inferenceWrapper(mcmcInfo)
    
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
%         tic
        waitbar(step/mcmcInfo.n_mcmc_steps,wb);     
        
        mcmcInfo.step = step;       
          
        % resample chains  promoter state sequences        
        mcmcInfo = resample_chains_v3(mcmcInfo);    % "Expectation Step"                 

        % perform additional "cross-talk" MH sampling if we are doing chain
        % tempering
%         if mcmcInfo.temperingFlag
%             mcmcInfo = temper_chains_v3(mcmcInfo);
%         end
                       
        % get empirical transition and occupancy counts    
        mcmcInfo = get_empirical_counts_v3(mcmcInfo);

        % use Gibbs sampling to update hyperparameters   
%         try
            mcmcInfo = update_hmm_parameters_v3(mcmcInfo);   
%         catch
%             break
%         end
        % if desired, perform MH moves to sample "nSteps"
        if mcmcInfo.inferNStepsFlag
            mcmcInfo = mh_sample_nSteps(mcmcInfo);           
%             if ~mcmcInfo.temperingFlag 
%                 mcmcInfo = mh_swap_nSteps(mcmcInfo);
%             end
        end
        
        % calculate updated logL
        mcmcInfo = calculateLogLikelihood(mcmcInfo);
%         toc
    end
    
    disp('done')
    delete(wb);