function log_fluo_diffs = calculate_fluo_logL(mcmcInfo)
                  
    % extrac useful parameters
    seq_length = mcmcInfo.seq_length;
    ind = mcmcInfo.ind;
    nSteps = mcmcInfo.nSteps;
    nStates = mcmcInfo.nStates;
    sigma = mcmcInfo.sigma;    
    coeff_MS2 = mcmcInfo.coeff_MS2;
    
    v_ref = repmat(reshape(mcmcInfo.v_curr,1,1,[]),1,mcmcInfo.n_chains);
    
    % get start and stop indices
    postInd = min([seq_length,ind+nSteps-1]);
    prevInd = max([1,ind-nSteps+1]);
    
    % extract relevant promoter state fragment
    initaition_fragment_init = mcmcInfo.v_curr(mcmcInfo.sample_chains_slice(prevInd:postInd,:));
    initiation_fragment = repmat(initaition_fragment_init,1,1,nStates);
    initiation_fragment(ind-prevInd+1,:,:) = v_ref;%repmat(mcmcInfo.sample_chains_slice(prevInd:postInd,:),1,1,nStates);           
                                        
    % calculate predicted fluorescence
    fluo_fragment = convn(coeff_MS2,initiation_fragment,'full');             
    fluo_fragment = fluo_fragment(ind-prevInd+1:end-length(coeff_MS2)+1,:,:);%+1(ind-prevInd+1:end-nSteps+1,:,:);
    
    % get differences
    log_fluo_diffs_full = ((mcmcInfo.observed_fluo(ind:postInd,mcmcInfo.trace_id)-fluo_fragment)./sigma).^2;
                                      
    % take average
    log_fluo_diffs = -mean(log_fluo_diffs_full,1);%-sum(0.5.*ms2_weights.*log_fluo_diffs_full,1)./ sum(coeff_MS2);% - log(sqrt(2*pi)*sigma);
    
    % reshape to be consistent with transition prob format    
    log_fluo_diffs = permute(log_fluo_diffs,[3 2 1]);
   