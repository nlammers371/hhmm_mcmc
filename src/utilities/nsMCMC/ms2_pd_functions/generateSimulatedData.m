function trueParams = generateSimulatedData(trueParams)
    
    % calculate MS2 convolution kernel
    trueParams.alpha = trueParams.alpha_frac*trueParams.nSteps;
    coeff_MS2 = ms2_loading_coeff_frac(trueParams.alpha, trueParams.nSteps, 12)';
    
    % calculate scaled noise
%     f_mean = (trueParams.pi0'*trueParams.v)*sum(coeff_MS2);   
%     trueParams.sigma = f_mean/10;        
    
    %%%%%%%%%%%%%%%%%%%%%%% simulate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trueParams.observed_fluo = NaN(trueParams.seq_length,trueParams.n_traces);
    trueParams.true_states = NaN(trueParams.seq_length,trueParams.n_traces);    
    trueParams.masterSimStruct = struct;
    for n = 1:trueParams.n_traces
        synthetic_data = synthetic_prob_frac(trueParams.seq_length, trueParams.alpha, trueParams.nStates, ...
                              trueParams.nSteps, trueParams.A, trueParams.v', trueParams.sigma, ...
                              trueParams.pi0);                                     

        trueParams.observed_fluo(:,n) = synthetic_data.fluo_MS2;
        trueParams.true_states(:,n) = synthetic_data.naive_states;
        % record full simulation info
        fieldNames = fieldnames(synthetic_data);
        for f = 1:length(fieldNames)
            trueParams.masterSimStruct(n).(fieldNames{f}) = synthetic_data.(fieldNames{f});
        end
    end       
    