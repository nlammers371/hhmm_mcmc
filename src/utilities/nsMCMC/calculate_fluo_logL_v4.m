function logL_fluo = calculate_fluo_logL_v4(mcmcInfo)
                  
    % note that this version uses normal array indexing
    
    % extract key hyperparameters
    nStates = mcmcInfo.nStates;
    n_chains = mcmcInfo.n_chains_eff;
    seq_length = mcmcInfo.seq_length;
    n_traces = mcmcInfo.n_traces;
    coeff_MS2 = mcmcInfo.MS2_kernel; 
    nStepsMax = size(coeff_MS2,1); 
    
    % calculate relevant indices 
    samp_index = mcmcInfo.samp_index;
    relevant_indices = mcmcInfo.step_ref + samp_index;
    n_pad_pre = sum(relevant_indices<=0);
    n_pad_post = sum(relevant_indices>seq_length);
    relevant_indices = relevant_indices(relevant_indices>0&relevant_indices<=seq_length);
    ind_filter = relevant_indices == samp_index;
    comp_filter = relevant_indices>=samp_index;
    comp_indices = relevant_indices(comp_filter);
    
    % extract useful parameters    
    sample_chains_curr = mcmcInfo.sample_chains_temp(2:end-1,:,:); % removing dummy start and ending states added for transition prob inference
    
    % generate array with state options 
    sample_chain_array = repmat(sample_chains_curr(relevant_indices,:,:),1,1,1,nStates);
    sample_chain_array(ind_filter,:,:,:) = repmat(reshape(1:nStates,1,1,1,nStates),1,n_chains,n_traces);   
    
    % generate array with corresponding initiation rates
    v_lin_indices = (sample_chain_array-1)*n_chains + mcmcInfo.chain_id_ref + 1;
    initiation_fragment = mcmcInfo.v_curr(v_lin_indices);
     
    % pad with zeros if necessary
%     initiation_fragment = initiation_fragment;
    if n_pad_pre > 0
        initiation_fragment = cat(1,zeros(n_pad_pre,n_chains,n_traces,nStates),initiation_fragment);
    end
%     if n_pad_post > 0
%         initiation_fragment = cat(1,initiation_fragment,zeros(n_pad_post,n_chains,n_traces,nStates));
%     end
    % use this to generate predicted fluorescence
%     fluo_fragment = NaN(size(initiation_fragment,1) + size(coeff_MS2,1)-1,n_chains,n_traces,nStates);
%     for i = 1:n_chains
%         fluo_fragment(:,i,:,:) = convn(coeff_MS2(:,i),initiation_fragment(:,i,:,:),'full');
%     end
    fluo_fragment = zeros(nStepsMax-n_pad_post,n_chains,n_traces,nStates);
    for i = 1:nStepsMax-n_pad_post
        fluo_fragment(i,:,:,:) = sum(coeff_MS2.*initiation_fragment(i:i+nStepsMax-1,:,:,:),1);
    end
    % cut back down to size
%     fluo_fragment = fluo_fragment(1:end-size(coeff_MS2,1)+1,:,:,:);
%     fluo_fragment = fluo_fragment(comp_filter,:,:,:);
    
    % extract corresponding fragment from experimental traces   
    if ~mcmcInfo.bootstrapFlag
        ref_fluo = repmat(permute(mcmcInfo.observed_fluo(comp_indices,:),[1 3 2]),1,n_chains,1,nStates);
    else
        ref_fluo = repmat(mcmcInfo.observed_fluo(comp_indices,:,:),1,1,1,nStates);
    end
       
    % generate sigma array    
    try
        sigma_ref = repmat(mcmcInfo.sigma_curr',size(ref_fluo,1),1,n_traces,nStates);
        logL_fluo_full = -0.5*(((ref_fluo-fluo_fragment)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));
    catch
        error('sigh')
    end
    % take average    
    logL_fluo = sum(logL_fluo_full,1);%-sum(0.5.*ms2_weights.*log_fluo_diffs_full,1)./ sum(coeff_MS2);% - log(sqrt(2*pi)*sigma);
   
    % reshape to be consistent with transition prob format    
    logL_fluo = permute(logL_fluo,[4 2 3 1]);
   