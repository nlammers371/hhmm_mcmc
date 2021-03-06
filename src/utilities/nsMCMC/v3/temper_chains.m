function mcmcInfo = temper_chains(mcmcInfo)

% extract parameters
A_log = log(mcmcInfo.A_curr(:,:,mcmcInfo.refChainVec));
coeff_MS2 = mcmcInfo.coeff_MS2;
pi0 = mcmcInfo.pi0_curr;
nStates = mcmcInfo.nStates;
n_traces = mcmcInfo.n_traces;
n_chains = mcmcInfo.n_chains;
nSteps = mcmcInfo.nSteps;
seq_length = mcmcInfo.seq_length;

% set options
n_rs_per_trace = mcmcInfo.n_rs_per_trace;

for c = mcmcInfo.n_chains-1:-1:1
    for t = 1:mcmcInfo.n_traces
      
        ref_trace = mcmcInfo.observed_fluo(:,t);            
        % extract chains
        t1 = mcmcInfo.sample_chains(:,c,t);
        t2 = mcmcInfo.sample_chains(:,c+1,t);
        
        % randomly propose swap(s)
        swap_options = find(t2~=t1);%2:length(t1)-1;
        swap_options = swap_options(~ismember(swap_options,[1 length(t1)]));
        n_samp = min([n_rs_per_trace,length(swap_options)]);
        if n_samp>1
            prop_points = randsample(swap_options,n_samp);
        else
            prop_points = swap_options;
        end
            
        % test each point
        for p = 1:length(prop_points)                    
            % make swap proposals
            t1_prop = t1;
            t1_prop(prop_points(p)) = t2(prop_points(p));
            t2_prop = t2;
            t2_prop(prop_points(p)) = t1(prop_points(p));
            
            % generate temporary array for indexing
            chain_array_temp = [t1 t2 t1_prop t2_prop];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Transition component 
            
            %%% previous state %%%
            % calculate linear indices 
            prev_state_array = chain_array_temp(prop_points(p)-1,:);
            curr_state_array = chain_array_temp(prop_points(p),:);
            row_col_array_from =  (prev_state_array-1)*nStates+curr_state_array;    

            % extract probabilities
            prev_probs_log = A_log(row_col_array_from);

            %%% following state %%%
            % calculate linear indices 
            post_state_array = chain_array_temp(prop_points(p)+1,:);
            row_col_array_to = (curr_state_array-1)*nStates+post_state_array;

            % extract probabilities
            post_probs_log = A_log(row_col_array_to);

            % combine
            logL_tr = prev_probs_log + post_probs_log;
          
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Emission component 
            %%%%%%%%%%%%%%%%%%%%%%%%%                        
            linear_indices = (chain_array_temp-1)*n_chains + [c c+1 c c+1];
            initiation_rates = mcmcInfo.v_curr(linear_indices);    

            fluo_predicted = convn(coeff_MS2,initiation_rates,'full');             
            fluo_predicted = fluo_predicted(1:end-length(coeff_MS2)+1,:,:);
                        
            sigma_ref = mcmcInfo.sigma_curr([c c+1 c c+1]);    
            
            logL_fluo = -0.5*(((ref_trace-fluo_predicted)./sigma_ref).^2 + log(2*pi*sigma_ref.^2));
                      
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Combine
            %%%%%%%%%%%%%%%%%%%%%%%%%             
            total_log_likelihoods = logL_tr + sum(logL_fluo);            
%             total_log_likelihoods = total_log_likelihoods./ mcmcInfo.tempGradVec([c c+1 c c+1]);
            temp1 = mcmcInfo.tempGradVec(c);
            temp2 = mcmcInfo.tempGradVec(c+1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % MH move            
            %%%%%%%%%%%%%%%%%%%%%%%%% 
            L_factor = exp((total_log_likelihoods(3)-total_log_likelihoods(1))/temp1)*...
                               exp((total_log_likelihoods(4)-total_log_likelihoods(2))/temp2);
%             move_flag = min([1,]);
            if L_factor > rand() 
                mcmcInfo.move_flag_array(p,c,t,mcmcInfo.step) = 1;
                mcmcInfo.sample_chains(:,c,t) = t1_prop;
                mcmcInfo.sample_chains(:,c+1,t) = t2_prop;
            end
        end
    end
end