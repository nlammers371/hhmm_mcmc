function mcmcInfo = predict_fluo_full_v3(mcmcInfo)
    
    coeff_MS2 = mcmcInfo.coeff_MS2_us;     
    n_chains = mcmcInfo.n_chains;
        
    linear_indices = (mcmcInfo.sample_chains-1)*n_chains + (1:n_chains);
    initiation_rates = mcmcInfo.v_curr(linear_indices)/mcmcInfo.upsample_factor;    

    fluo_out = fluo_conv_fun(initiation_rates,coeff_MS2);
    
    mcmcInfo.sample_fluo = fluo_out;