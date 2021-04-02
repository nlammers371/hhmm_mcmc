function mcmcInfo = genericInitialization(mcmcInfo)
    
    % calculate MS2 convolution kernel
    mcmcInfo.coeff_MS2 = ms2_loading_coeff(mcmcInfo.alpha, mcmcInfo.nSteps)';
    
    % calculate scaled noise
    f_mean = (mcmcInfo.trueParams.pi0'*mcmcInfo.trueParams.v)*sum(mcmcInfo.coeff_MS2);   
    mcmcInfo.trueParams.sigma = f_mean/10;
    
    %%%%%%%%%%%%%%%%%%%%%%% Generate helper arrays %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    mcmcInfo.state_options = 1:mcmcInfo.nStates;
    mcmcInfo.state_ref = repmat(reshape(mcmcInfo.state_options,1,1,[]),1,mcmcInfo.n_chains);

    % initialize arrays to store inference results
    n_updates = ceil(mcmcInfo.n_mcmc_steps/mcmcInfo.update_increment);
    n_chains = mcmcInfo.n_chains;
    mcmcInfo.logL_vec = NaN(n_updates,n_chains);
    mcmcInfo.A_inf_array = NaN(mcmcInfo.nStates,mcmcInfo.nStates,n_updates,n_chains);
    mcmcInfo.v_inf_array = NaN(n_updates,mcmcInfo.nStates,n_chains);
    mcmcInfo.pi0_inf_array = NaN(n_updates,mcmcInfo.nStates,n_chains);
    mcmcInfo.sigma_inf_array = NaN(n_updates,n_chains);

    %%%%%%%%%%%%%%%%%%%%%%% simulate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mcmcInfo.observed_fluo =cell(1,mcmcInfo.n_traces);
    trueParams = mcmcInfo.trueParams;
    mcmcInfo.masterSimStruct = struct;
    for n = 1:mcmcInfo.n_traces
        synthetic_data = synthetic_prob(mcmcInfo.seq_length, mcmcInfo.alpha, mcmcInfo.nStates, ...
                              mcmcInfo.nSteps, trueParams.A, trueParams.v', trueParams.sigma, trueParams.pi0);                                     

        mcmcInfo.observed_fluo{n} = synthetic_data.fluo_MS2';
        mcmcInfo.seq_len_vec(n) = length(mcmcInfo.observed_fluo{n});
        % record full simulation info
        fieldNames = fieldnames(synthetic_data);
        for f = 1:length(fieldNames)
            mcmcInfo.masterSimStruct(n).(fieldNames{f}) = synthetic_data.(fieldNames{f});
        end
    end       
        
    % generate reference vector
    mcmcInfo.chain_id_ref = 0:n_chains-1;
    mcmcInfo.temp_id_ref = reshape(0:mcmcInfo.n_temps_per_chain-1,1,1,[]);
    mcmcInfo.row_ref = (1:mcmcInfo.nStates)';
    mcmcInfo.step_ref = (-mcmcInfo.nSteps+1:mcmcInfo.nSteps-1)';
    
    % start parallel pool
    pool = gcp('nocreate');
    numcores = feature('numcores');
    targetCoreNum = min([numcores mcmcInfo.n_traces]);
    if isempty(pool) && targetCoreNum>4 %NL: this does not appear to be worthwhile for fewer than ~5 CPUs 
      parpool(targetCoreNum);      
    elseif ~isempty(pool) &&targetCorNum > 4
      delete(pool);
      parpool(targetCoreNum);            
    end