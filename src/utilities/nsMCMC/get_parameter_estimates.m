function [true_val_vec, param_mean_vec, param_ste_vec, param_names] = ...
                get_parameter_estimates(trueParams,A_inf_array,v_inf_array,rm_outlier_flag,inference_error_flag)
              
    % get number of states          
    nStates = trueParams.nStates;          
    n_chains = size(v_inf_array,3);
    % reorder states to be in ascending order as necessary
    for c = 1:n_chains
        v_slice = v_inf_array(:,:,c);
        v_mean = mean(v_slice,1);
        [~,si] = sort(v_mean);
        v_inf_array(:,:,c) = v_slice(:,si);
        
        a_slice = A_inf_array(:,:,:,c);
        A_inf_array(:,:,:,c) = a_slice(si,si,:);
    end
    
    % flag outlier chains
    outlier_chain_flags = inference_error_flag;
    if rm_outlier_flag
        % look for outliers in 2 key transition probs
        a21_vec = reshape(nanmean(A_inf_array(2,1,:,:),3),1,[]);
        [~,a21_flags] = rmoutliers(a21_vec);        
        a12_vec = reshape(nanmean(A_inf_array(1,2,:,:),3),1,[]);
        [~,a12_flags] = rmoutliers(a12_vec);
        
        % look for outliers in nonzero initiation rates
        v2_vec = reshape(nanmean(v_inf_array(:,2,:),1),1,[]);
        [~,v2_flags] = rmoutliers(v2_vec);
        v3_flags = false(size(v2_flags));
        if nStates == 3
            v3_vec = reshape(nanmean(v_inf_array(:,3,:),1),1,[]);
            [~,v3_flags] = rmoutliers(v3_vec);
        end
        outlier_chain_flags = outlier_chain_flags | a21_flags | a12_flags | v2_flags | v3_flags;        
    end
    
    % get transition parameter averages
    a_mean_vec = [];
    a_true_vec = [];
    a_ste_vec = [];
    a_string_cell = {};
    for i = 1:nStates
        j_vec = find(~ismember(1:nStates,i));
        for j = j_vec
              a_chunk = A_inf_array(j,i,:,~outlier_chain_flags);              
              
              a_mean_vec(end+1) = nanmean(a_chunk(:));
              a_ste_vec(end+1) = nanstd(a_chunk(:));
              
              a_string_cell(end+1) = {['a' num2str(j) num2str(i)]};
              a_true_vec(end+1) = trueParams.A(j,i);
        end
    end
    
    % get initiation parameter averages
    v_mean_vec = [];
    v_ste_vec = [];
    v_true_vec = [];
    v_string_cell = {};
    for i = 2:nStates        
          v_chunk = v_inf_array(:,i,~outlier_chain_flags);

          v_mean_vec(end+1) = nanmean(v_chunk(:));
          v_ste_vec(end+1) = nanstd(v_chunk(:));         
          
          v_string_cell(end+1) = {['v' num2str(i)]};        
          v_true_vec(end+1) = trueParams.v(i);
    end
    
    % combine
    param_mean_vec = [a_mean_vec v_mean_vec];
    param_ste_vec = [a_ste_vec v_ste_vec];
    param_names = [a_string_cell v_string_cell];
    true_val_vec = [a_true_vec v_true_vec];