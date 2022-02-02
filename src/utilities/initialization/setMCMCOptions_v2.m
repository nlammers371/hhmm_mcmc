function mcmcInfo = setMCMCOptions_v2(mcmcInfo, n_chains, inferMemory)

    mcmcInfo.n_chains = n_chains; % number of parallel MCMC chains to run   
    mcmcInfo.update_increment = 1;        
    
    %%%%%%%%%%%%%%%%%%%%% MCMC parameters %%%%%%%%%%%%%%%%  

    % memory parameter
    mcmcInfo.inferNStepsFlag = inferMemory;
    mcmcInfo.nStepsPropSize = 0.15;    
    mcmcInfo.nStep_tries_per_run = 10;
    mcmcInfo.nStepsGuess = 5 + rand()*3;
    mcmcInfo.nStepsMax = 10; % set upper limit
    mcmcInfo.nStepsMin = 3.5; % lower limit   

    % inference type
    mcmcInfo.refChainVec = true(1,mcmcInfo.n_chains_eff);                    
    mcmcInfo.move_flag_array = false(mcmcInfo.n_rs_per_trace,floor(n_chains*n_temps/2),mcmcInfo.n_traces,mcmcInfo.n_mcmc_steps);

    
        