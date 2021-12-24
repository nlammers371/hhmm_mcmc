function mcmcInfo = initializeInferenceArrays(mcmcInfo)

    %%%%%%%%%%%%%%%%%%%%%%% Generate helper arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    mcmcInfo.state_options = 1:mcmcInfo.nStates;
    mcmcInfo.state_ref = repmat(reshape(mcmcInfo.state_options,1,1,[]),1,mcmcInfo.n_chains);    
    
    % initialize arrays to store inference results
    n_updates = ceil(mcmcInfo.n_mcmc_steps/mcmcInfo.update_increment);
    n_chains = mcmcInfo.n_chains_eff;
    mcmcInfo.logL_vec = NaN(n_updates,n_chains);
    mcmcInfo.A_inf_array = NaN(mcmcInfo.nStates,mcmcInfo.nStates,n_updates,n_chains);
    mcmcInfo.v_inf_array = NaN(n_updates,mcmcInfo.nStates,n_chains);
    mcmcInfo.pi0_inf_array = NaN(n_updates,mcmcInfo.nStates,n_chains);
    mcmcInfo.sigma_inf_array = NaN(n_updates,n_chains);
    if mcmcInfo.inferNStepsFlag
        mcmcInfo.n_steps_inf_array = NaN(n_updates,n_chains);
    end
    
    % A prior--assume strongly diagonal PDF given short timescale
    % take A columns to follow multinomial Dirichlet distribution
    mcmcInfo.A_curr = NaN(mcmcInfo.nStates,mcmcInfo.nStates,mcmcInfo.n_chains*mcmcInfo.n_temps_per_chain);
    mcmcInfo.pi0_curr = NaN(mcmcInfo.n_chains*mcmcInfo.n_temps_per_chain,mcmcInfo.nStates);
    mcmcInfo.A_alpha = ones(mcmcInfo.nStates,mcmcInfo.nStates,mcmcInfo.n_chains*mcmcInfo.n_temps_per_chain);%*n_particles*n_traces;