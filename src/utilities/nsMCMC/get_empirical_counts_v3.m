function mcmcInfo = get_empirical_counts_v3(mcmcInfo)
    
    % extract parameters
    A_curr = mcmcInfo.A_curr;
    nStates = size(A_curr,1);
    n_chains = mcmcInfo.n_chains;

    % initialize count arrays 
    mcmcInfo.transition_count_array = zeros(size(A_curr,1),size(A_curr,2),n_chains);    
    mcmcInfo.initial_state_array = zeros(nStates,n_chains);   
    
    %%% record transitions
    
    % get linear indices
    from_array = mcmcInfo.sample_chains(1:end-1,:,:);
    to_array = mcmcInfo.sample_chains(2:end,:,:);        
    row_col_array = (from_array-1)*nStates + to_array;% + nStates^2*mcmcInfo.chain_id_ref;
    lin_index_array = row_col_array;
    
    % get state counts
    for k = 1:nStates
        mcmcInfo.initial_state_array(k,:) = sum(mcmcInfo.sample_chains(1,:,:)==k,3);% / n_chains;
    end        
        
    % get transition counts
    unique_indices = unique(lin_index_array(:));    
  
    n_vec = (0:n_chains-1)*nStates^2;
    for i = 1:length(unique_indices)        
        mcmcInfo.transition_count_array(n_vec+unique_indices(i)) = sum(sum(lin_index_array==unique_indices(i),1),3);
    end    
    
    mcmcInfo.state_counts = sum(mcmcInfo.transition_count_array,1);
%     if mcmcInfo.rateSamplingFlag&&mcmcInfo.adjustSamplingFlag % calculate adjusted counts if necessary
%         mcmcInfo = adjust_tr_counts(mcmcInfo);
%     end
    