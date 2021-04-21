function mcmcInfo = predict_fluo_full_v3(mcmcInfo)

    n_traces = mcmcInfo.n_traces;    
    coeff_MS2 = mcmcInfo.coeff_MS2;
    n_chains = mcmcInfo.n_chains_eff;
        
    linear_indices = (mcmcInfo.sample_chains-1)*n_chains + (1:n_chains);
    initiation_rates = mcmcInfo.v_curr(linear_indices);    

    fluo_fragment = NaN(size(initiation_rates,1)+size(coeff_MS2,1)-1,n_chains,n_traces);
    for n = 1:n_chains
        fluo_fragment(:,n,:) = convn(coeff_MS2(:,n),initiation_rates(:,n,:),'full');             
    end
    mcmcInfo.sample_fluo = fluo_fragment(1:end-size(coeff_MS2,1)+1,:,:);