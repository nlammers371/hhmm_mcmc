function [logL_fluo, new_fluo_array, comp_indices] = calculate_fluo_logL_v5(mcmcInfo)
                      
    
    % extract key hyperparameters
    nStates = mcmcInfo.nStates;
    n_chains = mcmcInfo.n_chains_eff;   
    n_traces = mcmcInfo.n_traces;
%     coeff_MS2 = mcmcInfo.MS2_kernel; 
    us_factor = mcmcInfo.upsample_factor;
    seq_length = mcmcInfo.seq_length*us_factor;
    coeff_MS2_us = mcmcInfo.MS2_kernel_us; 
%     if ~mcmcInfo.inferNStepsFlag
%         coeff_MS2 = coeff_MS2(:,1);
%     end
    nStepsMax = size(coeff_MS2_us,1);         
    
    % calculate relevant indices 
    samp_index = mcmcInfo.samp_index;
    relevant_indices =  (samp_index:samp_index+nStepsMax-1)';
    comp_filter = relevant_indices<=seq_length;
    comp_indices = relevant_indices(comp_filter);    
    comp_indices_fluo = ceil((comp_indices(1))/us_factor):floor((comp_indices(end))/us_factor); % map to experimental data indices
    coeff_MS2 = coeff_MS2_us(comp_filter,:);
    us_map_filter = ismember(comp_indices/us_factor,comp_indices_fluo);
    v_curr_long = mcmcInfo.v_curr_long(comp_filter,:,:);
    
    % extract w X n_chains X n_traces piece of fluo array
    fluo_array = mcmcInfo.sample_fluo_temp(comp_indices,:,:);
    
    % extract 1 X n_chains X n_traces array of current states
    current_states = mcmcInfo.sample_chains_temp(samp_index+1,:,:);
    
    % get corresponding emission rate array
    v_lin_indices = mcmcInfo.chain_id_ref + 1 + (current_states-1)*n_chains;
    initiation_fragment = mcmcInfo.v_curr(v_lin_indices)/us_factor;

    % calculate expected contribution of current state to each time step
    fluo_curr = initiation_fragment.*coeff_MS2;        
    
    % extract relevant comparison fragment
    ref_fluo = permute(mcmcInfo.observed_fluo(comp_indices_fluo,:,:),[1 3 2]);
    
    % generate sigma array       
    sigma_ref = repmat(mcmcInfo.sigma_curr',size(ref_fluo,1),1,n_traces);
    
    % now iterate through each possible state and calculate expected
    % change, new predicted fluo, and resulting likelihood
    logL_fluo = NaN(nStates,n_chains,n_traces);
    new_fluo_array = NaN(size(coeff_MS2,1),n_chains,n_traces,nStates);
    
    for k = 1:nStates
        fluo_alt = v_curr_long(:,:,k).*coeff_MS2/us_factor;        
        fluo_new = fluo_array + fluo_alt - fluo_curr;
        fluo_cp = fluo_new(us_map_filter,:,:);
        
        % calculate likelihood        
        logL_fluo_full = -0.5*(((ref_fluo-fluo_cp)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));
        logL_fluo(k,:,:) = sum(logL_fluo_full,1);
        
        % store fluo
        new_fluo_array(:,:,:,k) = fluo_new;
    end    
   