function log_fluo_diffs = calculate_fluo_logL_full(mcmcInfo)
                  
    % extrac useful parameters
    seq_length = mcmcInfo.seq_length;
    ind = mcmcInfo.ind;
    nSteps = mcmcInfo.nSteps;
    nStates = mcmcInfo.nStates;
    sigma = mcmcInfo.sigma;
    state_ref = mcmcInfo.state_ref;
    coeff_MS2 = mcmcInfo.coeff_MS2;
    
    v_ref = repmat(reshape(mcmcInfo.v_curr,1,1,[]),1,mcmcInfo.n_chains);
    
    % get start and stop indices
    postInd = min([seq_length,ind+nSteps-1]);
    prevInd = max([1,ind-nSteps+1]);
    
    % extract relevant promoter state fragment
    initaition_fragment_init = repmat(mcmcInfo.v_curr(mcmcInfo.sample_chains_slice),1,1,nStates);%repmat(mcmcInfo.sample_chains_slice(prevInd:postInd,:),1,1,nStates);   
    initaition_fragment_init(ind,:,:) = v_ref;
    initiation_fragment = initaition_fragment_init;%(prevInd:postInd,:,:);%vertcat(initaition_fragment_init(prevInd:ind-1,:,:), v_ref,...
                           %       initaition_fragment_init(ind+1:postInd,:,:));
                                    
    
    % calculate predicted fluorescence
    fluo_fragment = convn(coeff_MS2,initiation_fragment,'full');             
    fluo_fragment = fluo_fragment(1:end-length(coeff_MS2)+1,:,:);%+1(ind-prevInd+1:end-nSteps+1,:,:);
    
    % get differences
    log_fluo_diffs_full = ((mcmcInfo.true_fluo-...(ind:postInd)-...
                                        fluo_fragment)./sigma).^2;
                                      
    % take weighted average
    ms2_weights = coeff_MS2(1:postInd-ind+1);
    
    log_fluo_diffs = -mean(log_fluo_diffs_full(ind:postInd,:,:),1);%-sum(0.5.*ms2_weights.*log_fluo_diffs_full,1)./ sum(coeff_MS2);% - log(sqrt(2*pi)*sigma);
    
    % reshape to be consistent with transition prob format
    try
        log_fluo_diffs = reshape(log_fluo_diffs, nStates,[]);
    catch
        error('wtf')
    end