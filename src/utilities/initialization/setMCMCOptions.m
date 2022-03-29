function mcmcInfo = setMCMCOptions(mcmcInfo, trueParams)
    
    % Set defaults
    defaultOptions = struct;
    defaultOptions.ensembleInferenceFlag = 0;
    defaultOptions.strongAPriorFlag = 0;
    defaultOptions.bootstrapFlag = 0;
    defaultOptions.bootstrapFlagPar = 0;
    defaultOptions.annealingSigmaFlag = 0;
    defaultOptions.ensembleInferenceFlag = 0;
    defaultOptions.mhResamplingFlag = 0;
    defaultOptions.inferNStepsFlag = 0;
    defaultOptions.n_chains = 25;
    defaultOptions.update_increment = 1; 
    defaultOptions.upsample_factor = 1; 
    fnamesDefault = fieldnames(defaultOptions);
    
    % check to see if any of the above options have been specified
    fnamesInput = fieldnames(mcmcInfo);
    optionFlags = ~ismember(fnamesDefault,fnamesInput);
    for i = find(optionFlags)'
        mcmcInfo.(fnamesDefault{i}) = defaultOptions.(fnamesDefault{i});
    end                   
    
    % multiple jumps per time point is currently only supported with MH
    % algorithm
    if mcmcInfo.upsample_factor > 1
        mcmcInfo.mhResamplingFlag = 1;
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
        mcmcInfo.observed_fluo = NaN(mcmcInfo.seq_length_data, mcmcInfo.n_chains, mcmcInfo.n_traces);
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

    mcmcInfo.n_chains_eff = mcmcInfo.n_chains;

    % memory parameter    
    mcmcInfo.nStepsPropSize = 0.15;    
    mcmcInfo.nStep_tries_per_run = 10;
    mcmcInfo.nStepsGuess = (5 + rand()*3) * trueParams.upsample_factor;
    mcmcInfo.nStepsMax = 12 * trueParams.upsample_factor; % set upper limit
    mcmcInfo.nStepsMin = 3.5 * trueParams.upsample_factor; % lower limit   

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
    mcmcInfo.chainIDVec = 1:mcmcInfo.n_chains;
    
    % transfer fields
    mcmcInfo.nStates = trueParams.nStates;    
    mcmcInfo.alpha_frac = trueParams.alpha_frac;
    mcmcInfo.observed_fluo = trueParams.observed_fluo;
    if ~defaultOptions.inferNStepsFlag
        mcmcInfo.nSteps_true = trueParams.nSteps_true;
        mcmcInfo.nSteps_data = trueParams.nSteps_data;
    end    


    
        