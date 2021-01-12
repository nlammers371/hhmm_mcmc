function mcmcInfo = predict_fluo_full_ens(mcmcInfo)

%     n_traces = mcmcInfo.n_traces;    
    coeff_MS2 = mcmcInfo.coeff_MS2;
%     n_chains = mcmcInfo.n_chains;
    
    % convert to emission values
%     emissions_full = mcmcInfo.v_curr(mcmcInfo.sample_chains);
%     mcmcInfo.sample_fluo = NaN(size(mcmcInfo.sample_chains));
        
    initiation_rates = mcmcInfo.v_curr(mcmcInfo.sample_chains);    
    
%     for n = 1:n_traces
%         emissions = emissions_full(:,:,n);                
%         
%         % convolve with emissions vector
%         fluo_new = conv2(emissions,coeff_MS2,'full');
%         mcmcInfo.sample_fluo(:,:,n) = fluo_new(1:end-length(coeff_MS2)+1,:);
%         
%     end
    
    fluo_fragment = convn(coeff_MS2,initiation_rates,'full');             
    mcmcInfo.sample_fluo = fluo_fragment(1:end-length(coeff_MS2)+1,:,:);