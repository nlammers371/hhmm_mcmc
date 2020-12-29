function mcmcInfo = predict_fluo_full(mcmcInfo)

    n_traces = mcmcInfo.n_traces;    
    coeff_MS2 = mcmcInfo.coeff_MS2;    
    n_chains = mcmcInfo.n_chains;
    seq_length = mcmcInfo.seq_length;
    v_curr = mcmcInfo.v_curr;
    
    % convert to emission values
    if ~mcmcInfo.par_chain_flag
        initiation_events_full = v_curr(mcmcInfo.sample_chains);
    else        
        row_index_vec = repmat(repelem(1:n_chains,1,seq_length)',n_traces,1);        
        linear_indices = sub2ind(size(v_curr),row_index_vec,mcmcInfo.sample_chains(:));
        initiation_events_full = reshape(v_curr(linear_indices),seq_length,n_chains,n_traces);
    end        
    % convolve with emissions vector
    fluo_new = convn(coeff_MS2,initiation_events_full ,'full'); 

    mcmcInfo.sample_fluo = fluo_new(1:end-length(coeff_MS2)+1,:,:);            