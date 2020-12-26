function mcmcInfo = predict_fluo_full(mcmcInfo)

    n_traces = mcmcInfo.n_traces;    
    coeff_MS2 = mcmcInfo.coeff_MS2;
    
    % convert to emission values
    emissions_full = mcmcInfo.v_curr(mcmcInfo.sample_chains);
    mcmcInfo.sample_fluo = NaN(size(mcmcInfo.sample_chains));
    
    for n = 1:n_traces
        emissions = emissions_full(:,:,n);                
        
        % convolve with emissions vector
        fluo_new = conv2(emissions,coeff_MS2,'full');
        mcmcInfo.sample_fluo(:,:,n) = fluo_new(1:end-length(coeff_MS2)+1,:);
        
    end