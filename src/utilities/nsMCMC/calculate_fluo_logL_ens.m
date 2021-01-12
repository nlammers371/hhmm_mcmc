function logL_fluo = calculate_fluo_logL_ens(mcmcInfo)
                  
    % extrac useful parameters
    seq_length_dummy = size(mcmcInfo.sample_chains_dummy,1);%mcmcInfo.seq_length;
    seq_length = mcmcInfo.seq_length;
    indexArray = mcmcInfo.indArray;
    nSteps = mcmcInfo.nSteps;
    nStates = mcmcInfo.nStates;
    n_traces = mcmcInfo.n_traces;
    n_chains = mcmcInfo.n_chains;
    sigma_curr = mcmcInfo.sigma_curr;    
    coeff_MS2 = mcmcInfo.coeff_MS2;
    
%     v_ref = permute(repmat(mcmcInfo.v_curr',1,n_chains,n_traces,1),[4,2,3,1]);
    
    % get start and stop indices    
    indexArrayFull = indexArray + mcmcInfo.step_ref + nSteps - 1;
    linIndexArrayFull = indexArrayFull + mcmcInfo.chain_id_ref*seq_length_dummy + mcmcInfo.trace_id_ref*seq_length_dummy*n_chains;  
    
    % extract states
    initiation_fragment_init = mcmcInfo.v_curr(mcmcInfo.sample_chains_dummy(linIndexArrayFull));
    
    % zero out dummy placeholder states
    dummyFilter = indexArrayFull-nSteps+1<1 | indexArrayFull-nSteps+1>seq_length;
    initiation_fragment_init(dummyFilter) = 0;
    
    % add proposed states    
    initiation_fragment = repmat(initiation_fragment_init,1,1,1,nStates);
    initiation_fragment(nSteps,:,:,:) = repmat(reshape(mcmcInfo.v_curr,1,1,1,[]),1,n_chains,n_traces);%repmat(mcmcInfo.sample_chains_slice(prevInd:postInd,:),1,1,nStates);           
                                        
    % calculate predicted fluorescence
    fluo_fragment = convn(coeff_MS2,initiation_fragment,'full');             
    fluo_fragment = fluo_fragment(nSteps:2*nSteps-1,:,:,:);
    
    % get differences   
    linIndexArrayFluo = indexArrayFull(nSteps:end,:,:) +  mcmcInfo.trace_id_ref*seq_length_dummy;
    ref_fluo = repmat(mcmcInfo.observed_fluo_dummy(linIndexArrayFluo),1,1,1,nStates);
    ref_fluo(dummyFilter(nSteps:end,:,:)) = fluo_fragment(dummyFilter(nSteps:end,:,:));
    % set leading and trailing observations equal to model fluo. This is
    % the easiest way to remove from likelihood calculation
    
    logL_fluo_full = -0.5*(((ref_fluo-fluo_fragment)./sigma_curr).^2 + log(2*pi*sigma_curr.^2));
                                      
    % take average    
    logL_fluo = sum(logL_fluo_full,1);%-sum(0.5.*ms2_weights.*log_fluo_diffs_full,1)./ sum(coeff_MS2);% - log(sqrt(2*pi)*sigma);
   
    % reshape to be consistent with transition prob format    
    logL_fluo = permute(logL_fluo,[4 2 3 1]);
   