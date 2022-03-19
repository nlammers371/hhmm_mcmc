function mcmcInfo = setMCMCOptions(mcmcInfo, n_chains, inferMemory)

    mcmcInfo.n_chains = n_chains; % number of parallel MCMC chains to run   
    mcmcInfo.update_increment = 1;        
    
    %%%%%%%%%%%%%%%%%%%%% MCMC parameters %%%%%%%%%%%%%%%%
%     if ~temperingFlag
%         n_temps = 1;
%     end
    
    % tempering options
%     mcmcInfo.temperingFlag = false;%temperingFlag; % use parallel tempering
%     mcmcInfo.n_rs_per_trace = n_swaps*temperingFlag; % number of swap proposals per neighboring trace pair
%     mcmcInfo.n_temps_per_chain = n_temps; % number of rungs in the temperature ladder for each chain
    mcmcInfo.n_chains_eff = mcmcInfo.n_chains;%mcmcInfo.n_temps_per_chain*mcmcInfo.n_chains;

    % memory parameter
    mcmcInfo.inferNStepsFlag = inferMemory;
    mcmcInfo.nStepsPropSize = 0.15;    
    mcmcInfo.nStep_tries_per_run = 10;
    mcmcInfo.nStepsMax = 12*mcmcInfo.upsample_factor; % set upper limit
    mcmcInfo.nStepsMin = 3.5*mcmcInfo.upsample_factor; % lower limit   

    % set curve if we are using sigma-mediated annealing 
    if mcmcInfo.annealingSigmaFlag
        stepVec = 1:mcmcInfo.n_mcmc_steps;
        sg = mcmcInfo.burn_in;
        gradientVec = 2*exp(-stepVec/sg*5)+1;
        mcmcInfo.tempGradArray = repmat(gradientVec',1,mcmcInfo.n_chains);
        mcmcInfo.tempGradArray(mcmcInfo.burn_in+1:end,:) = 1;
    end
    
    % inference type
    mcmcInfo.refChainVec = true(1,mcmcInfo.n_chains_eff);
%     mcmcInfo.refChainVec(1:mcmcInfo.n_temps_per_chain:end) = true; % designate T=1 chains that will actually be used for inference
    mcmcInfo.chainIDVec = 1:mcmcInfo.n_chains;%mcmcInfo.n_temps_per_chain);
%     mcmcInfo.infoSharingFlag = info_sharing;
%     mcmcInfo.temp_incrememt = temp_increment;%sqrt(2); % from: https://emcee.readthedocs.io/en/v2.2.1/user/pt/
%     exp_vec = repmat(0:mcmcInfo.n_temps_per_chain-1,1,mcmcInfo.n_chains);
%     mcmcInfo.tempGradVec = mcmcInfo.temp_incrememt.^exp_vec;%logspace(0,log10(mcmcInfo.max_temp),mcmcInfo.n_chains);

    % mcmcInfo.tempGradVec(2) = 1;
%     mcmcInfo.move_flag_array = false(mcmcInfo.n_rs_per_trace,floor(n_chains*n_temps/2),mcmcInfo.n_traces,mcmcInfo.n_mcmc_steps);

    
        