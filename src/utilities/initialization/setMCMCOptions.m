function mcmcInfo = setMCMCOptions(mcmcInfo)
    
    % set defaults
    defaultOptions = struct;
%     defaultOptions.ensembleInferenceFlag = 0;
    defaultOptions.strongAPriorFlag = 0;
%     defaultOptions.bootstrapFlag = 0;
%     defaultOptions.bootstrapFlagPar = 0;
    defaultOptions.annealingSigmaFlag = 0;
%     defaultOptions.ensembleInferenceFlag = 0;
    defaultOptions.mhResamplingFlag = 0;
%     defaultOptions.PFResamplingFlag = 0;
    defaultOptions.upsample_factor = 1;
    
    defaultOptions.rateSamplingFlag = 0;
    defaultOptions.adjustSamplingFlag = 0;
    defaultOptions.rateSamplingHRFlag = 0;
    defaultOptions.mhQSamplingFlag = 0;
%     defaultOptions.mhInferenceFlag = 0;
%     defaultOptions.reducedModelFlag = 0;
    
    defaultOptions.inferNStepsFlag = 0;
    defaultOptions.n_chains = 25;
    defaultOptions.update_increment = 1; 
    defaultOptions.save_trace_results = 0; 
    defaultOptions.em_timer_flag = 0;
    fnamesDefault = fieldnames(defaultOptions);
    
    % check to see if any of the above options have been specified
    fnamesInput = fieldnames(mcmcInfo);
    optionFlags = ~ismember(fnamesDefault,fnamesInput);
    for i = find(optionFlags)'
        mcmcInfo.(fnamesDefault{i}) = defaultOptions.(fnamesDefault{i});
    end               
        
    %%%%%%%%%%%%%%%%%%%%% MCMC parameters %%%%%%%%%%%%%%%%
    mcmcInfo.trace_id_array = repmat((1:mcmcInfo.n_traces)',1,mcmcInfo.n_chains);   
    mcmcInfo.n_traces_per_chain = mcmcInfo.n_traces;    

    % memory parameter
    if mcmcInfo.inferNStepsFlag
        mcmcInfo.nStepsPropSize = 0.15;    
        mcmcInfo.nStep_tries_per_run = 10;
    end
    mcmcInfo.nStepsGuess = 5 + rand()*3;
    mcmcInfo.nStepsMax = 20; % set upper limit
    mcmcInfo.nStepsMin = 3; % lower limit   5

    % Rate matrix sampling options
    if mcmcInfo.mhQSamplingFlag
        mcmcInfo.QPropSize = 0.005;    
        mcmcInfo.Q_n_tries_per_run = 1;        
    end
    mcmcInfo.QMax = 0.2;
    
    % set curve if we are using sigma-mediated annealing 
    if mcmcInfo.annealingSigmaFlag
        stepVec = 1:mcmcInfo.n_mcmc_steps;
        sg = mcmcInfo.burn_in;
        gradientVec = 2*exp(-stepVec/sg*5)+1;
        mcmcInfo.tempGradArray = repmat(gradientVec',1,mcmcInfo.n_chains);
        mcmcInfo.tempGradArray(mcmcInfo.burn_in+1:end,:) = 1;
    end
    
    % inference type
    mcmcInfo.refChainVec = true(1,mcmcInfo.n_chains);
    mcmcInfo.chainIDVec = 1:mcmcInfo.n_chains;

    
        