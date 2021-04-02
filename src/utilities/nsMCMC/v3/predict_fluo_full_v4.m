function mcmcInfo = predict_fluo_full_v4(mcmcInfo)

%     n_traces = mcmcInfo.n_traces;    
    coeff_MS2 = mcmcInfo.coeff_MS2;
    n_chains = mcmcInfo.n_chains;
    n_traces = mcmcInfo.n_traces;   
    
    % convert to emission values
%     emissions_full = mcmcInfo.v_curr(mcmcInfo.sample_chains);
%     mcmcInfo.sample_fluo = NaN(size(mcmcInfo.sample_chains));
    for m = 1:n_traces
        linear_indices = (mcmcInfo.sample_chains{m}-1)*n_chains + (1:n_chains);
        initiation_rates = mcmcInfo.v_curr(linear_indices);    

        fluo_fragment = convn(coeff_MS2,initiation_rates,'full');             
        mcmcInfo.sample_fluo{m} = fluo_fragment(1:end-length(coeff_MS2)+1,:,:);
    end