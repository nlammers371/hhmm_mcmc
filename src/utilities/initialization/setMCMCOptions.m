function mcmcInfo = setMCMCOptions(mcmcInfo)
    
    % set defaults
    defaultOptions = struct;
    defaultOptions.ensembleInferenceFlag = 0;
    defaultOptions.strongAPriorFlag = 0;
    defaultOptions.bootstrapFlag = 0;
    defaultOptions.bootstrapFlagPar = 0;
    defaultOptions.annealingSigmaFlag = 0;
    defaultOptions.ensembleInferenceFlag = 0;
    defaultOptions.mhResamplingFlag = 0;
    defaultOptions.PFResamplingFlag = 0;
    defaultOptions.upsample_factor = 1;
    defaultOptions.rateSamplingFlag = 0;
    
    defaultOptions.mhInferenceFlag = 0;
    defaultOptions.reducedModelFlag = 0;
    
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

    if mcmcInfo.bootstrapFlag
        % randomly assign subset of traces to each chain
%         mcmcInfo.n_traces_per_chain = ceil(200/mcmcInfo.seq_length);
        mcmcInfo.trace_id_array = NaN(mcmcInfo.n_traces_per_chain,mcmcInfo.n_chains);
        for n = 1:mcmcInfo.n_chains
            mcmcInfo.trace_id_array(:,n) = randsample(1:mcmcInfo.n_traces,mcmcInfo.n_traces_per_chain,true);
        end
        
        % now update key variables to be compatible with downstream
        % functions
        mcmcInfo.n_traces_orig = mcmcInfo.n_traces;
        mcmcInfo.n_traces = mcmcInfo.n_traces_per_chain;
        mcmcInfo.observed_fluo_orig = mcmcInfo.observed_fluo;
        mcmcInfo.observed_fluo = NaN(mcmcInfo.seq_length, mcmcInfo.n_chains, mcmcInfo.n_traces);
        for n = 1:n_chains
            mcmcInfo.observed_fluo(:,n,:) = mcmcInfo.observed_fluo_orig(:,mcmcInfo.trace_id_array(:,n)');
        end
    elseif mcmcInfo.bootstrapFlagPar        
        mcmcInfo.trace_id_array = randsample(1:mcmcInfo.n_traces,mcmcInfo.n_traces_per_chain,true); 
        mcmcInfo.n_traces_orig = mcmcInfo.n_traces;
        mcmcInfo.n_traces = mcmcInfo.n_traces_per_chain;     
        mcmcInfo.observed_fluo_orig = mcmcInfo.observed_fluo;
        mcmcInfo.observed_fluo = mcmcInfo.observed_fluo_orig(:,mcmcInfo.trace_id_array');
        
    else
        mcmcInfo.n_traces_per_chain = mcmcInfo.n_traces;
    end
%     if ~temperingFlag
%         n_temps = 1;
%     end
    
    % tempering options
%     mcmcInfo.temperingFlag = false;%temperingFlag; % use parallel tempering
%     mcmcInfo.n_rs_per_trace = n_swaps*temperingFlag; % number of swap proposals per neighboring trace pair
%     mcmcInfo.n_temps_per_chain = n_temps; % number of rungs in the temperature ladder for each chain
    mcmcInfo.n_chains_eff = mcmcInfo.n_chains;%mcmcInfo.n_temps_per_chain*mcmcInfo.n_chains;

    % memory parameter
    mcmcInfo.nStepsPropSize = 0.15;    
    mcmcInfo.nStep_tries_per_run = 10;
    mcmcInfo.nStepsGuess = 5 + rand()*3;
    mcmcInfo.nStepsMax = 20; % set upper limit
    mcmcInfo.nStepsMin = 3; % lower limit   
%     if ~mcmcInfo.inferNStepsFlag 
%         mcmcInfo.nStepsMax = ceil(mcmcInfo.nSteps);
%     end
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

    
        