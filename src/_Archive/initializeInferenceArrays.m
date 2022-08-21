function mcmcInfo = initializeInferenceArrays(mcmcInfo)

    %%%%%%%%%%%%%%%%%%%%%%% Generate helper arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    mcmcInfo.state_options = 1:mcmcInfo.nStates;
    mcmcInfo.state_ref = repmat(reshape(mcmcInfo.state_options,1,1,[]),1,mcmcInfo.n_chains);    
    
    % initialize arrays to store inference results
    n_updates = ceil(mcmcInfo.n_mcmc_steps/mcmcInfo.update_increment);
    n_chains = mcmcInfo.n_chains_eff;
    mcmcInfo.logL_vec = NaN(n_updates,n_chains);
    
    if ~mcmcInfo.reducedModelFlag && ~mcmcInfo.mhInferenceFlag
        mcmcInfo.A_inf_array = NaN(mcmcInfo.nStates,mcmcInfo.nStates,n_updates,n_chains);
        mcmcInfo.v_inf_array = NaN(n_updates,mcmcInfo.nStates,n_chains);
        mcmcInfo.pi0_inf_array = NaN(n_updates,mcmcInfo.nStates,n_chains);
        mcmcInfo.sigma_inf_array = NaN(n_updates,n_chains);
        if mcmcInfo.inferNStepsFlag
            mcmcInfo.n_steps_inf_array = NaN(n_updates,n_chains);
        end
    elseif mcmcInfo.mhInferenceFlag && mcmcInfo.reducedModelFlag
        mcmcInfo.k_inf_array = NaN(n_updates,mcmcInfo.nStates,1);
        mcmcInfo.r_inf_array = NaN(n_updates,mcmcInfo.nStates-1,1);
        mcmcInfo.sigma_inf_array = NaN(n_updates,1);
        if mcmcInfo.inferNStepsFlag
            mcmcInfo.n_steps_inf_array = NaN(n_updates,1);
        end
    else
        error('Incompatible inference options')
    end
           
    % initialize arrays to store current parameters
    mcmcInfo.A_curr = NaN(mcmcInfo.nStates,mcmcInfo.nStates,mcmcInfo.n_chains);
    mcmcInfo.pi0_curr = NaN(mcmcInfo.n_chains_eff,mcmcInfo.nStates);
    mcmcInfo.v_curr = NaN(mcmcInfo.n_chains_eff,mcmcInfo.nStates);   
    mcmcInfo.sigma_curr = NaN(mcmcInfo.n_chains,1);
    
    if mcmcInfo.em_timer_flag
        mcmcInfo.em_time_vec = NaN(mcmcInfo.n_mcmc_steps,1);
    end
   
      