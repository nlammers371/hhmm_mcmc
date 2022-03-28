function [true_val_vec, param_mean_vec, param_ste_vec, param_names] = ...
                get_parameter_estimates_boot(trueParams,A_inf_array,v_inf_array,logL_vec,boot_id_vec,inference_error_flags)
              
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
    
    % identify most likely chain for each bootstrap
    ml_id_vec = [];
    boot_index = unique(boot_id_vec);
    for b = 1:length(boot_index)
        boot_indices = find(boot_id_vec==boot_index(b)&~inference_error_flags);
        [~,ml_i] = nanmax(logL_vec(boot_indices));
        ml_id_vec(end+1) = boot_indices(ml_i);
    end
    
    % subset to only ml chains
    A_inf_array_ml = A_inf_array(:,:,:,ml_id_vec);
    v_inf_array_ml = v_inf_array(:,:,ml_id_vec);
    
    % get transition parameter averages
    a_mean_vec = [];
    a_true_vec = [];
    a_ste_vec = [];
    a_string_cell = {};
    for i = 1:nStates
        j_vec = find(~ismember(1:nStates,i));
        for j = j_vec
              a_chunk = A_inf_array_ml(j,i,:,:);              
              
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
          v_chunk = v_inf_array_ml(:,i,:);

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