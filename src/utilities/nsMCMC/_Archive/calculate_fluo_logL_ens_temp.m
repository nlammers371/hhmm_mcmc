function logL_fluo = calculate_fluo_logL_ens_temp(mcmcInfo)
                  
    % extrac useful parameters
    seq_length = mcmcInfo.seq_length;
    ind = mcmcInfo.ind;
    nSteps = mcmcInfo.nSteps;
    nStates = mcmcInfo.nStates;
    n_traces = mcmcInfo.n_traces;
    n_chains = mcmcInfo.n_chains;
    sigma_curr = repelem(mcmcInfo.sigma_curr,n_chains)';    
    coeff_MS2 = mcmcInfo.coeff_MS2;
    
    v_curr = repmat(mcmcInfo.v_curr',n_chains,1);
    v_ref = permute(repmat(v_curr',1,1,n_traces,1),[4,2,3,1]);
    
    % get start and stop indices
    postInd = min([seq_length,ind+nSteps-1]);
    prevInd = max([1,ind-nSteps+1]);
    
    % extract relevant promoter state fragment   mcmcInfo.masterSimStruct.naive_states'
    linear_indices = (mcmcInfo.sample_chains(prevInd:postInd,:,:)-1)*n_chains + (1:n_chains);
    initiation_fragment_init = v_curr(linear_indices);
    initiation_fragment = repmat(initiation_fragment_init,1,1,1,nStates);

    initiation_fragment(ind-prevInd+1,:,:,:) = v_ref;%repmat(mcmcInfo.sample_chains_slice(prevInd:postInd,:),1,1,nStates);           

    % calculate predicted fluorescence
    fluo_fragment = convn(coeff_MS2,initiation_fragment,'full');             
    fluo_fragment = fluo_fragment(ind-prevInd+1:end-length(coeff_MS2)+1,:,:,:);
    
    % get differences   
    ref_fluo = repmat(permute(mcmcInfo.observed_fluo(ind:postInd,:,:),[1 3 2]),1,n_chains,1,nStates);
    sigma_ref = repmat(sigma_curr',size(ref_fluo,1),1,n_traces,nStates);
    
    logL_fluo_full = -0.5*(((ref_fluo-fluo_fragment)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));
                                      
    % take average    
    logL_fluo = sum(logL_fluo_full,1);%-sum(0.5.*ms2_weights.*log_fluo_diffs_full,1)./ sum(coeff_MS2);% - log(sqrt(2*pi)*sigma);
   
    % reshape to be consistent with transition prob format    
    logL_fluo = permute(logL_fluo,[4 2 3 1]);
   