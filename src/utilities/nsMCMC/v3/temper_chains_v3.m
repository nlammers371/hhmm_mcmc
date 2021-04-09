function mcmcInfo = temper_chains_v3(mcmcInfo)

% extract parameters
A_log = log(mcmcInfo.A_curr);
coeff_MS2 = mcmcInfo.coeff_MS2;

% pi0 = mcmcInfo.pi0_curr;
nStates = mcmcInfo.nStates;
n_traces = mcmcInfo.n_traces;
n_chains = mcmcInfo.n_chains;

% nSteps = mcmcInfo.nSteps;
seq_length = mcmcInfo.seq_length;
n_temps_per_chain = mcmcInfo.n_temps_per_chain; 

% set options
n_rs_per_trace = mcmcInfo.n_rs_per_trace;
n_chains_eff = mcmcInfo.n_chains_eff; 
chain_id_eff_ref = 0:n_chains_eff-1;

%% generate reference arrays

id_vec = 1:n_temps_per_chain*n_chains;
temp_id_vec = repmat(1:n_temps_per_chain,1,n_chains);
chain_id_vec = repelem(1:n_chains,1,n_temps_per_chain);
n_swaps = floor(n_chains*n_temps_per_chain/2);
index_vec = 1:seq_length;

% iterate through traces
for m = 1:n_traces
  
    % extract relevant slice
    chain_slice = mcmcInfo.sample_chains(:,:,m);

    % extract trace
    ref_trace = mcmcInfo.observed_fluo(:,m);   
    
    % iterate through swaps
    for rs = 1:n_rs_per_trace
        %%%%%%%%%%%%%%%
        % randomly assign non-ovoerlaping swap pairs
        %%%%%%%%%%%%%%%
        
        % initialize arrays
        trace_array = NaN(seq_length,n_swaps,2);
        temperature_array = NaN(2,n_swaps);
        id_array = NaN(2,n_swaps);  
        chain_id_array = NaN(1,n_swaps,2);  
        used_id_flags = false(size(id_vec));
        
        % randomly assign partners
        for n = 1:n_swaps          
          
            % pick first
            id1 = randsample(find(~used_id_flags),1);
            used_id_flags(id1) = 1;

            % find permitted partners
            chain_neighbor_options = chain_id_vec == chain_id_vec(id1) & abs(temp_id_vec-temp_id_vec(id1))==1; 
            temp_neighbor_options = temp_id_vec == temp_id_vec(id1) & abs(chain_id_vec-chain_id_vec(id1))==1; 
            partner_options = (temp_neighbor_options | chain_neighbor_options) & ~used_id_flags; 
            if any(partner_options) 
                partner_ids = repelem(find(partner_options),2);% keep ransample from being stupid            
            else
                partner_ids = repelem(find(~used_id_flags),2);
            end

            % select second id
            id2 = randsample(partner_ids,1);
            used_id_flags(id2) = 1;

            % assign
            id_array(:,n) = [id1 id2];
            trace_array(:,n,:) = chain_slice(:,[id1 id2]);
            temperature_array(:,n) = mcmcInfo.tempGradVec([id1 id2]);
            chain_id_array(1,n,:) = chain_id_eff_ref([id1 id2]);
        end
        
        % propose time points to swap
        mismatch_flags = 1*(trace_array(:,:,1)~=trace_array(:,:,2)) + 1e-10;
        
        % exclude start and end points for now
        mismatch_flags(1,:) = 0;
        mismatch_flags(end,:) = 0;
        
        % sample
        swap_points = NaN(1,size(mismatch_flags,2));
        for i = 1:size(mismatch_flags,2)
            swap_points(i) = randsample(index_vec,1,true,mismatch_flags(:,i));
        end
        
        %%%%%%%%%%%%%%%%%
        % perform MH swap move 
        
        % linearize
        swap_lin_indices = (0:n_swaps-1)*seq_length + swap_points;
        swap_ind_array = swap_lin_indices + reshape((0:3)*seq_length*n_swaps,1,1,[]);
        
        % create temp array
        trace_array_prop = trace_array;
        trace_array_prop(swap_lin_indices) = trace_array(swap_lin_indices+seq_length*n_swaps);
        trace_array_prop(swap_lin_indices+seq_length*n_swaps) = trace_array(swap_lin_indices);
        
        trace_array_full = cat(3,trace_array,trace_array_prop);
        chain_id_array_full = repmat(chain_id_array,1,1,2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Transition component 

        %%% previous state %%%
        % calculate linear indices 
        prev_state_array = trace_array_full(swap_ind_array-1);
        curr_state_array = trace_array_full(swap_ind_array);
        row_col_array_from = chain_id_array_full.*nStates^2 + (prev_state_array-1)*nStates + curr_state_array;    

        % extract probabilities
        prev_probs_log = A_log(row_col_array_from);

        %%% following state %%%
        % calculate linear indices 
        post_state_array = trace_array_full(swap_ind_array+1);
        row_col_array_to = chain_id_array_full.*nStates^2 + (curr_state_array-1)*nStates + post_state_array;

        % extract probabilities
        post_probs_log = A_log(row_col_array_to);

        % combine
        logL_tr = prev_probs_log + post_probs_log;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Emission component 
        %%%%%%%%%%%%%%%%%%%%%%%%%                        
        linear_indices = (trace_array_full-1)*n_chains_eff + chain_id_array_full+1;
        initiation_rates = mcmcInfo.v_curr(linear_indices);    

        fluo_predicted = convn(coeff_MS2,initiation_rates,'full');             
        fluo_predicted = fluo_predicted(1:end-length(coeff_MS2)+1,:,:);

        sigma_ref = mcmcInfo.sigma_curr(chain_id_array_full+1);    

        logL_fluo = -0.5*(((ref_trace-fluo_predicted)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));
        
        
        % combine
        total_log_likelihoods = logL_tr + sum(logL_fluo);  
        
        % calculate MH metric
        L_factor = exp((total_log_likelihoods(1,:,3)-total_log_likelihoods(1,:,1))./temperature_array(1,:) + ...
                        (total_log_likelihoods(1,:,4)-total_log_likelihoods(1,:,2))./temperature_array(2,:));
                      
        % perform MH move
        mh_flags = L_factor >= rand(size(L_factor));
        
        % update
        trace_array(swap_lin_indices(mh_flags)) = trace_array_prop(swap_lin_indices(mh_flags));
        trace_array(swap_lin_indices(mh_flags)+seq_length*n_swaps) = trace_array_prop(swap_lin_indices(mh_flags)+seq_length*n_swaps);
        
        % map back to main array
        swap_ids = find(mh_flags);
        for i = swap_ids
            chain_slice(swap_points(i),id_array(:,i)) = trace_array(swap_points(i),i,:);
        end
    end
    mcmcInfo.sample_chains(:,:,m) = chain_slice;
    
end



% %%
% for c = mcmcInfo.n_chains-1:-1:1
%     for t = 1:mcmcInfo.n_traces
%       
%         ref_trace = mcmcInfo.observed_fluo(:,t);            
%         % extract chains
%         t1 = mcmcInfo.sample_chains(:,c,t);
%         t2 = mcmcInfo.sample_chains(:,c+1,t);
%         
%         % randomly propose swap(s)
%         swap_options = find(t2~=t1);%2:length(t1)-1;
%         swap_options = swap_options(~ismember(swap_options,[1 length(t1)]));
%         n_samp = min([n_rs_per_trace,length(swap_options)]);
%         if n_samp>1
%             prop_points = randsample(swap_options,n_samp);
%         else
%             prop_points = swap_options;
%         end
%             
%         % test each point
%         for p = 1:length(prop_points)                    
%             % make swap proposals
%             t1_prop = t1;
%             t1_prop(prop_points(p)) = t2(prop_points(p));
%             t2_prop = t2;
%             t2_prop(prop_points(p)) = t1(prop_points(p));
%             
%             % generate temporary array for indexing
%             chain_array_temp = [t1 t2 t1_prop t2_prop];
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%
%             % Transition component 
%             
%             %%% previous state %%%
%             % calculate linear indices 
%             prev_state_array = chain_array_temp(prop_points(p)-1,:);
%             curr_state_array = chain_array_temp(prop_points(p),:);
%             row_col_array_from =  (prev_state_array-1)*nStates+curr_state_array + chain_id_array;    
% 
%             % extract probabilities
%             prev_probs_log = A_log(row_col_array_from);
% 
%             %%% following state %%%
%             % calculate linear indices 
%             post_state_array = chain_array_temp(prop_points(p)+1,:);
%             row_col_array_to = (curr_state_array-1)*nStates+post_state_array;
% 
%             % extract probabilities
%             post_probs_log = A_log(row_col_array_to);
% 
%             % combine
%             logL_tr = prev_probs_log + post_probs_log;
%           
%             %%%%%%%%%%%%%%%%%%%%%%%%%
%             % Emission component 
%             %%%%%%%%%%%%%%%%%%%%%%%%%                        
%             linear_indices = (chain_array_temp-1)*n_chains + [c c+1 c c+1];
%             initiation_rates = mcmcInfo.v_curr(linear_indices);    
% 
%             fluo_predicted = convn(coeff_MS2,initiation_rates,'full');             
%             fluo_predicted = fluo_predicted(1:end-length(coeff_MS2)+1,:,:);
%                         
%             sigma_ref = mcmcInfo.sigma_curr([c c+1 c c+1]);    
%             
%             logL_fluo = -0.5*(((ref_trace-fluo_predicted)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));
%                       
%             %%%%%%%%%%%%%%%%%%%%%%%%%
%             % Combine
%             %%%%%%%%%%%%%%%%%%%%%%%%%             
%             total_log_likelihoods = logL_tr + sum(logL_fluo);            
% %             total_log_likelihoods = total_log_likelihoods./ mcmcInfo.tempGradVec([c c+1 c c+1]);
%             temp1 = mcmcInfo.tempGradVec(c);
%             temp2 = mcmcInfo.tempGradVec(c+1);
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%
%             % MH move            
%             %%%%%%%%%%%%%%%%%%%%%%%%% 
%             L_factor = exp((total_log_likelihoods(3)-total_log_likelihoods(1))/temp1)*...
%                                exp((total_log_likelihoods(4)-total_log_likelihoods(2))/temp2);
% %             move_flag = min([1,]);
%             if L_factor > rand() 
%                 mcmcInfo.move_flag_array(p,c,t,mcmcInfo.step) = 1;
%                 mcmcInfo.sample_chains(:,c,t) = t1_prop;
%                 mcmcInfo.sample_chains(:,c+1,t) = t2_prop;
%             end
%         end
%     end
% end