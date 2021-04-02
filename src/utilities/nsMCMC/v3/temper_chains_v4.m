function mcmcInfo = temper_chains_v4(mcmcInfo)

% extract parameters
A_log = log(mcmcInfo.A_curr);
coeff_MS2 = mcmcInfo.coeff_MS2;
pi0 = mcmcInfo.pi0_curr;
nStates = mcmcInfo.nStates;
n_traces = mcmcInfo.n_traces;
n_chains = mcmcInfo.n_chains;
chain_id_ref = mcmcInfo.chain_id_ref;
nSteps = mcmcInfo.nSteps;
% seq_length = mcmcInfo.seq_length;

% set options
n_rs_per_trace = mcmcInfo.n_rs_per_trace;

% propagate info down temperature ladder
for temp = mcmcInfo.n_temps_per_chain-1:-1:1
    for m = 1:mcmcInfo.n_traces
         
        % extract trace
        ref_trace = mcmcInfo.observed_fluo{m}; 
        seq_len = length(ref_trace);
        
        % extract chains
        temp1_trace = mcmcInfo.sample_chains{m}(:,:,temp);
        temp2_trace = mcmcInfo.sample_chains{m}(:,:,temp+1);
        
        % randomly propose swap(s)
        swap_options = 2:size(temp1_trace,1)-1; % exclude endpoints for now
        prop_points_array = reshape(randsample(swap_options,n_rs_per_trace*n_chains),n_rs_per_trace,n_chains);
        
        % test each point
        for p = 1:n_rs_per_trace    
          
            % generate linear indices
            prop_points = prop_points_array(p,:);
            prop_points_lin = prop_points + chain_id_ref*seq_len;
            
            % make swap proposals
            temp1_trace_prop = temp1_trace;
            temp1_trace_prop(prop_points_lin) = temp2_trace(prop_points_lin);
            temp2_trace_prop = temp2_trace;
            temp2_trace_prop(prop_points_lin) = temp1_trace(prop_points_lin);
            
            % generate temporary array for indexing
            chain_array_temp = cat(3,temp1_trace, temp2_trace, temp1_trace_prop, temp2_trace_prop);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Transition component of logL
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%% previous state %%%
            % calculate linear indices 
            prop_points_lin_3D = prop_points_lin + reshape((0:3)*numel(temp1_trace),1,1,4);
            
            prev_state_array = chain_array_temp(prop_points_lin_3D-1);
            curr_state_array = chain_array_temp(prop_points_lin_3D);
            row_col_array_from =  (prev_state_array-1)*nStates+curr_state_array;    

            % extract probabilities
            prev_probs_log = A_log(row_col_array_from);

            %%% following state %%%
            % calculate linear indices 
            post_state_array = chain_array_temp(prop_points_lin_3D+1);
            row_col_array_to = (curr_state_array-1)*nStates+post_state_array;

            % extract probabilities
            post_probs_log = A_log(row_col_array_to);

            % combine
            logL_tr = prev_probs_log + post_probs_log;
          
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Emission logL component 
            %%%%%%%%%%%%%%%%%%%%%%%%%     
            
            linear_indices = (chain_array_temp-1)*n_chains + repmat((mcmcInfo.chain_id_ref)+1,1,1,4);%[temp temp+1 temp temp+1];
            initiation_rates = mcmcInfo.v_curr(linear_indices);    

            fluo_predicted = convn(coeff_MS2,initiation_rates,'full');             
            fluo_predicted = fluo_predicted(1:end-length(coeff_MS2)+1,:,:);
                        
            sigma_ref = mcmcInfo.sigma_curr;    
            
            logL_fluo = -0.5*(((ref_trace-fluo_predicted)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Combine
            %%%%%%%%%%%%%%%%%%%%%%%%%             
            total_log_likelihoods = logL_tr + sum(logL_fluo);            
            total_log_likelihoods = total_log_likelihoods./ reshape(mcmcInfo.tempGradVec([temp temp+1 temp temp+1]),1,1,4);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % MH move            
            %%%%%%%%%%%%%%%%%%%%%%%%% 
            logL_factor = exp(sum(total_log_likelihoods(:,:,3:4),3)-sum(total_log_likelihoods(:,:,1:2),3));
%             move_flag = min([1,]);
            move_flags = logL_factor > rand(size(logL_factor));
            % update where appropriate
            sample_slice1 = mcmcInfo.sample_chains{m}(:,:,temp);
            sample_slice2 = mcmcInfo.sample_chains{m}(:,:,temp+1);
            % update
            sample_slice1(prop_points_lin(move_flags)) = temp1_trace_prop(prop_points_lin(move_flags));
            sample_slice2(prop_points_lin(move_flags)) = temp2_trace_prop(prop_points_lin(move_flags));
            
            mcmcInfo.sample_chains{m}(:,:,temp) = sample_slice1;
            mcmcInfo.sample_chains{m}(:,:,temp+1) = sample_slice2;
%             end
        end
    end
end