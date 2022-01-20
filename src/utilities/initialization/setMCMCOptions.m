function mcmcInfo = setMCMCOptions(mcmcInfo, n_chains, temperingFlag, n_temps, n_swaps, inferMemory, MCMCInitFlag, temp_increment, info_sharing)

    mcmcInfo.n_chains = n_chains; % number of parallel MCMC chains to run
   
    mcmcInfo.update_increment = 2;
    mcmcInfo.swapSteps = mcmcInfo.burn_in + mcmcInfo.swapInc + 1:mcmcInfo.swapInc:mcmcInfo.n_mcmc_steps;
    
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
    mcmcInfo.temperingFlag = temperingFlag; % use parallel tempering
    mcmcInfo.n_rs_per_trace = n_swaps*temperingFlag; % number of swap proposals per neighboring trace pair
    mcmcInfo.n_temps_per_chain = n_temps; % number of rungs in the temperature ladder for each chain
    mcmcInfo.n_chains_eff = mcmcInfo.n_temps_per_chain*mcmcInfo.n_chains;
    mcmcInfo.enforceRateConsistency = 0;

    % memory parameter
    mcmcInfo.inferNStepsFlag = inferMemory;
    mcmcInfo.nStepsPropSize = 0.2;    
%     mcmcInfo.trueParams.nSteps = 6; % True parameters
    mcmcInfo.nStepsGuess = 4 + rand()*3;
    mcmcInfo.nStepsMax = 10; % set upper limit
    mcmcInfo.nStepsMin = 3.5; % lower limit   

    % inference type
    mcmcInfo.refChainVec = false(1,mcmcInfo.n_chains_eff);
    mcmcInfo.refChainVec(1:mcmcInfo.n_temps_per_chain:end) = true; % designate T=1 chains that will actually be used for inference
    mcmcInfo.chainIDVec = repelem(1:mcmcInfo.n_chains,mcmcInfo.n_temps_per_chain);
    mcmcInfo.infoSharingFlag = info_sharing;
    mcmcInfo.temp_incrememt = temp_increment;%sqrt(2); % from: https://emcee.readthedocs.io/en/v2.2.1/user/pt/
    exp_vec = repmat(0:mcmcInfo.n_temps_per_chain-1,1,mcmcInfo.n_chains);
    mcmcInfo.tempGradVec = mcmcInfo.temp_incrememt.^exp_vec;%logspace(0,log10(mcmcInfo.max_temp),mcmcInfo.n_chains);

    % mcmcInfo.tempGradVec(2) = 1;
    mcmcInfo.move_flag_array = false(mcmcInfo.n_rs_per_trace,floor(n_chains*n_temps/2),mcmcInfo.n_traces,mcmcInfo.n_mcmc_steps);

    
        