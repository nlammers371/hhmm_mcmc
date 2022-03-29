function trueParams = generateSimulatedData(trueParams)
    
    % calculate MS2 convolution kernel
    trueParams.alpha = trueParams.alpha_frac*trueParams.nSteps_true;    
    
    %%%%%%%%%%%%%%%%%%%%%%% simulate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trueParams.observed_fluo = NaN(trueParams.seq_length_true,trueParams.n_traces);
    trueParams.true_states = NaN(trueParams.seq_length_true,trueParams.n_traces);    
    trueParams.masterSimStruct = struct;
    for n = 1:trueParams.n_traces
        synthetic_data = synthetic_prob_frac(trueParams.seq_length_true, trueParams.alpha, trueParams.nStates, ...
                              trueParams.nSteps_true, trueParams.A, trueParams.v', trueParams.sigma, ...
                              trueParams.pi0);                                     

        trueParams.observed_fluo(:,n) = synthetic_data.fluo_MS2;
        trueParams.true_states(:,n) = synthetic_data.naive_states;
        % record full simulation info
        fieldNames = fieldnames(synthetic_data);
        for f = 1:length(fieldNames)
            trueParams.masterSimStruct(n).(fieldNames{f}) = synthetic_data.(fieldNames{f});
        end
    end       
    
    % downsample data    
    trueParams.observed_fluo_full = trueParams.observed_fluo;
    trueParams.observed_fluo = trueParams.observed_fluo(1:trueParams.upsample_factor:end,:);