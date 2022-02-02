function mcmcInfo = simulateRawData(mcmcInfo)             
    
    %%%%%%%%%%%%%%%%%%%%%%% simulate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mcmcInfo.observed_fluo = NaN(mcmcInfo.seq_length,mcmcInfo.n_traces);
    mcmcInfo.true_states = NaN(mcmcInfo.seq_length,mcmcInfo.n_traces);
    trueParams = mcmcInfo.trueParams;
    mcmcInfo.masterSimStruct = struct;
    for n = 1:mcmcInfo.n_traces
        synthetic_data = synthetic_prob_frac(mcmcInfo.seq_length, trueParams.alpha, mcmcInfo.nStates, ...
                              trueParams.nSteps, trueParams.A, trueParams.v', trueParams.sigma, ...
                              trueParams.pi0);                                     

        mcmcInfo.observed_fluo(:,n) = synthetic_data.fluo_MS2;
        mcmcInfo.true_states(:,n) = synthetic_data.naive_states;
        % record full simulation info
        fieldNames = fieldnames(synthetic_data);
        for f = 1:length(fieldNames)
            mcmcInfo.masterSimStruct(n).(fieldNames{f}) = synthetic_data.(fieldNames{f});
        end
    end       
    