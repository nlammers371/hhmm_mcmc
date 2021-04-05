function mcmcInfo = genericInitialization(mcmcInfo)
    
    % calculate MS2 convolution kernel
    mcmcInfo.coeff_MS2 = ms2_loading_coeff_frac(mcmcInfo.alpha, mcmcInfo.nStepsCurr, mcmcInfo.nStepsMax)';
    
    % calculate scaled noise
    f_mean = (mcmcInfo.trueParams.pi0'*mcmcInfo.trueParams.v)*sum(mcmcInfo.coeff_MS2);   
    mcmcInfo.trueParams.sigma = f_mean/10;
    
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

    %%%%%%%%%%%%%%%%%%%%%%% simulate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mcmcInfo.observed_fluo = NaN(mcmcInfo.seq_length,mcmcInfo.n_traces);
    trueParams = mcmcInfo.trueParams;
    mcmcInfo.masterSimStruct = struct;
    for n = 1:mcmcInfo.n_traces
        synthetic_data = synthetic_prob(mcmcInfo.seq_length, mcmcInfo.alpha, mcmcInfo.nStates, ...
                              trueParams.nSteps, trueParams.A, trueParams.v', trueParams.sigma, ...
                              trueParams.pi0);                                     

        mcmcInfo.observed_fluo(:,n) = synthetic_data.fluo_MS2;
        % record full simulation info
        fieldNames = fieldnames(synthetic_data);
        for f = 1:length(fieldNames)
            mcmcInfo.masterSimStruct(n).(fieldNames{f}) = synthetic_data.(fieldNames{f});
        end
    end       
    