function logL_fluo = calculate_fluo_logL_v3_beta(mcmcInfo)
                  
    % extrac useful parameters
    seq_length_dummy = size(mcmcInfo.sample_chains_dummy,1);%mcmcInfo.seq_length;
    seq_length_dummy2 = size(mcmcInfo.sample_fluo_dummy2,1);
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
    naive_states = mcmcInfo.sample_chains_dummy(linIndexArrayFull);
    v_lin_indices = (naive_states-1)*n_chains + mcmcInfo.chain_id_ref + 1;
    initiation_fragment_init = mcmcInfo.v_curr(v_lin_indices);
    
    % zero out dummy placeholder states
    dummyFilter = indexArrayFull-nSteps+1<1 | indexArrayFull-nSteps+1>seq_length;
    initiation_fragment_init(dummyFilter) = 0;
    
    % add proposed states    
    initiation_fragment = repmat(initiation_fragment_init,1,1,1,nStates);
    v_perm = permute(mcmcInfo.v_curr,[4,1,3,2]);
    initiation_fragment(nSteps,:,:,:) = repmat(v_perm,1,1,n_traces);%repmat(mcmcInfo.sample_chains_slice(prevInd:postInd,:),1,1,nStates);           
                                        
    % calculate predicted fluorescence
    fluo_fragment = convn(coeff_MS2,initiation_fragment,'full');             
    fluo_fragment = fluo_fragment(nSteps:2*nSteps-1,:,:,:);
    
    %%
    indexArrayTrunc = indexArray + mcmcInfo.step_ref_trunc;
    linIndexArrayTrunc = indexArrayTrunc + mcmcInfo.chain_id_ref*seq_length_dummy2 + mcmcInfo.trace_id_ref*seq_length_dummy2*n_chains;  
    predicted_fluo_curr = mcmcInfo.sample_fluo_dummy2(linIndexArrayTrunc);
    
    %%
    
    % get differences   
    linIndexArrayFluo = indexArrayFull(nSteps:end,:,:) +  mcmcInfo.trace_id_ref*seq_length_dummy;
    ref_fluo = repmat(mcmcInfo.observed_fluo_dummy(linIndexArrayFluo),1,1,1,nStates);
    dummy_trunc = dummyFilter(nSteps:end,:,:);
    % set leading and trailing observations equal to model fluo. This is
    % the easiest way to remove from likelihood calculation
    ref_fluo(dummy_trunc) = fluo_fragment(dummy_trunc);
       
    % generate sigma array
    sigma_ref = repmat(sigma_curr,nSteps,1,n_traces,nStates);
    logL_fluo_full = -0.5*(((ref_fluo-fluo_fragment)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));
 
    % take average    
    logL_fluo = sum(logL_fluo_full,1);%-sum(0.5.*ms2_weights.*log_fluo_diffs_full,1)./ sum(coeff_MS2);% - log(sqrt(2*pi)*sigma);
   
    % reshape to be consistent with transition prob format    
    logL_fluo = permute(logL_fluo,[4 2 3 1]);
   