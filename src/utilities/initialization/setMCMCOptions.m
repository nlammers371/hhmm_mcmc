function mcmcInfo = setMCMCOptions(mcmcInfo, n_chains, temperingFlag, n_temps, n_swaps, inferMemory, MCMCInitFlag)

    mcmcInfo.n_chains = n_chains; % number of parallel MCMC chains to run
   
    
    %%%%%%%%%%%%%%%%%%%%% MCMC parameters %%%%%%%%%%%%%%%%
    if ~temperingFlag
        n_temps = 1;
    end
    mcmcInfo.MCMCInitFlag = 0;
    if MCMCInitFlag
        mcmcInfo.MCMCInitFlag = 1;
        mcmcInfo.MCMCInitFactor = 10;
        mcmcInfo.n_chains = mcmcInfo.MCMCInitFactor*mcmcInfo.n_chains;
        mcmcInfo.n_mcmc_steps = 25;
        mcmcInfo.n_reps = 1;
    end
    
    % tempering options
    mcmcInfo.temperingFlag = temperingFlag; % use parallel tempering?
    mcmcInfo.n_rs_per_trace = n_swaps; % number of swap proposals per neighboring trace pair
    mcmcInfo.n_temps_per_chain = n_temps; % number of rungs in the temperature ladder for each chain
    mcmcInfo.n_chains_eff = mcmcInfo.n_temps_per_chain*mcmcInfo.n_chains;
    mcmcInfo.enforceRateConsistency = 0;

    % memory parameter
    mcmcInfo.inferNStepsFlag = inferMemory;
    mcmcInfo.nStepsPropSize = 0.1;    
%     mcmcInfo.trueParams.nSteps = 6; % True parameters
    mcmcInfo.nStepsGuess = 6;
    mcmcInfo.nStepsMax = 9; % set upper limit
    mcmcInfo.nStepsMin = 4; % lower limit   

    % inference type
    mcmcInfo.refChainVec = false(1,mcmcInfo.n_chains_eff);
    mcmcInfo.refChainVec(1:mcmcInfo.n_temps_per_chain:end) = true; % designate T=1 chains that will actually be used for inference
    mcmcInfo.chainIDVec = repelem(1:mcmcInfo.n_chains,mcmcInfo.n_temps_per_chain);
    mcmcInfo.temp_incrememt = sqrt(2); % from: https://emcee.readthedocs.io/en/v2.2.1/user/pt/
    exp_vec = repmat(0:mcmcInfo.n_temps_per_chain-1,1,mcmcInfo.n_chains);
    mcmcInfo.tempGradVec = mcmcInfo.temp_incrememt.^exp_vec;%logspace(0,log10(mcmcInfo.max_temp),mcmcInfo.n_chains);

    % mcmcInfo.tempGradVec(2) = 1;
    mcmcInfo.move_flag_array = false(mcmcInfo.n_rs_per_trace,mcmcInfo.n_chains-1,mcmcInfo.n_traces,mcmcInfo.n_mcmc_steps);

    
        