function trueParams = generateSimulatedData(trueParams)
    
    % calculate MS2 convolution kernel
    trueParams.alpha = trueParams.alpha_frac*trueParams.nSteps;
%     coeff_MS2 = ms2_loading_coeff_frac(trueParams.alpha, trueParams.nSteps, 12)';
    
    % calculate scaled noise
%     f_mean = (trueParams.pi0'*trueParams.v)*sum(coeff_MS2);   
%     trueParams.sigma = f_mean/10;        
    
    %%%%%%%%%%%%%%%%%%%%%%% simulate traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trueParams.observed_fluo = NaN(trueParams.seq_length,trueParams.n_traces);
    trueParams.true_states = NaN(trueParams.seq_length,trueParams.n_traces);    
    trueParams.masterSimStruct = struct;
    for n = 1:trueParams.n_traces
        if trueParams.discrete_data_flag 
            synthetic_data = synthetic_prob_frac(trueParams.seq_length, trueParams.alpha, trueParams.nStates, ...
                                  trueParams.nSteps, trueParams.A, trueParams.v', trueParams.sigma, ...
                                  trueParams.pi0);     
        else
            trueParams.r_emission = trueParams.v'/trueParams.tres;
            synthetic_data = synthetic_rate_gillespie_frac(trueParams.seq_length, trueParams.alpha, trueParams.nStates, ...
                                 trueParams.nSteps, trueParams.R, trueParams.tres, trueParams.r_emission, trueParams.sigma,...
                                 trueParams.pi0); 
        end

        trueParams.observed_fluo(:,n) = synthetic_data.fluo_MS2;
        if trueParams.discrete_data_flag 
            trueParams.true_states(:,n) = synthetic_data.naive_states;
        else
            trueParams.true_states(:,n) = interp1(synthetic_data.transition_times,synthetic_data.naive_states,0:trueParams.tres:trueParams.tres*(trueParams.seq_length-1),'previous');
        end
        % record full simulation info
        fieldNames = fieldnames(synthetic_data);
        for f = 1:length(fieldNames)
            trueParams.masterSimStruct(n).(fieldNames{f}) = synthetic_data.(fieldNames{f});
        end
    end       
    