function mcmcInfo = predict_fluo_full_v3(mcmcInfo)
    
    coeff_MS2 = mcmcInfo.coeff_MS2;     
    n_chains = mcmcInfo.n_chains_eff;
        
    linear_indices = (mcmcInfo.sample_chains-1)*n_chains + (1:n_chains);
    initiation_rates = mcmcInfo.v_curr(linear_indices);    

    fluo_out = fluo_conv_fun(initiation_rates,coeff_MS2);
    
    mcmcInfo.sample_fluo = fluo_out;