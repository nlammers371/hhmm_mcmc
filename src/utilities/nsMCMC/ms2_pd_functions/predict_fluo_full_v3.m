function mcmcInfo = predict_fluo_full_v3(mcmcInfo)
    
    coeff_MS2 = mcmcInfo.coeff_MS2;     
    n_chains = mcmcInfo.n_chains_eff;
        
    linear_indices = (mcmcInfo.sample_chains-1)*n_chains + (1:n_chains);
    initiation_rates = mcmcInfo.v_curr(linear_indices);    

    fluo_out = fluo_conv_fun(initiation_rates,coeff_MS2);
    
    if mcmcInfo.upsample_factor > 1
        mcmcInfo.sample_fluo_full = fluo_out;
        mcmcInfo.sample_fluo = fluo_out(1:mcmcInfo.upsample_factor:end,:,:);
    else
        mcmcInfo.sample_fluo = fluo_out;
    end