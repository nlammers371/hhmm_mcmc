function logL_fluo = calculate_fluo_logL_par(mcmcInfo)
                  
    % extrac useful parameters
    seq_length = mcmcInfo.seq_length;
    ind = mcmcInfo.ind;
    nSteps = mcmcInfo.nSteps;
    nStates = mcmcInfo.nStates;
    n_chains = mcmcInfo.n_chains;
    sigma_curr = mcmcInfo.sigma_curr;    
    coeff_MS2 = mcmcInfo.coeff_MS2;
    
    v_ref = mcmcInfo.v_curr';
    
    % get start and stop indices
    postInd = seq_length;%min([seq_length,ind+nSteps-1]);
    prevInd = 1;%max([1,ind-nSteps+1]);
    
    % extract relevant promoter state fragment    
    linear_indices = (mcmcInfo.sample_chains-1)*n_chains + (1:n_chains);
    initaition_fragment_init = mcmcInfo.v_curr(linear_indices);
    initiation_fragment = repmat(initaition_fragment_init,1,1,1,nStates);
    initiation_fragment(ind,:,:,:) = v_ref;%repmat(mcmcInfo.sample_chains_slice(prevInd:postInd,:),1,1,nStates);           
                                        
    % calculate predicted fluorescence
    fluo_fragment = convn(coeff_MS2,initiation_fragment,'full');             
    fluo_fragment = fluo_fragment(1:end-length(coeff_MS2)+1,:,:);
%     fluo_fragment = fluo_fragment(ind-prevInd+1:end-length(coeff_MS2)+1,:,:);%+1(ind-prevInd+1:end-nSteps+1,:,:);
    
    % get differences    
    logL_fluo_full = -0.5*((mcmcInfo.observed_fluo-fluo_fragment)./sigma_curr).^2 - log(sqrt(2*pi)*sigma_curr);
                                      
    % take average    
    logL_fluo = sum(logL_fluo_full,1);%-sum(0.5.*ms2_weights.*log_fluo_diffs_full,1)./ sum(coeff_MS2);% - log(sqrt(2*pi)*sigma);
   
    % reshape to be consistent with transition prob format    
    logL_fluo = permute(logL_fluo,[3 2 1]);
   