function mcmcInfo = update_hmm_parameters_v4(mcmcInfo)    

    % extrace parameters
    nStates = mcmcInfo.nStates;
    seq_length = mcmcInfo.seq_length;
    n_traces = mcmcInfo.n_traces;
    n_chains = mcmcInfo.n_chains;    
    coeff_MS2_us = mcmcInfo.coeff_MS2_us;        
    us_factor = mcmcInfo.upsample_factor; 
    
    %% %%%%%%%%%% update transition matrix (A) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    update_flag = mod(mcmcInfo.step,mcmcInfo.update_increment) == 0 || mcmcInfo.step == mcmcInfo.n_mcmc_steps;
    mcmcInfo.update_flag = update_flag;
    if mcmcInfo.update_increment~=1
        update_index = ceil(mcmcInfo.step/mcmcInfo.update_increment) + 1;
    else
        update_index = mcmcInfo.step;
    end
        
    if ~mcmcInfo.mhQSamplingFlag
        % perform Gibbs update
        mcmcInfo = gibbs_update_A(mcmcInfo);                
    else
        % perform MH update
        mcmcInfo = mh_update_A(mcmcInfo);
    end
                
    
    %% %%%%%%%%%%%%% update emission vector(V) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update V       
    if true%~mcmcInfo.adjustSamplingFlag
        [mcmcInfo, F_array, y_array] = get_state_counts(mcmcInfo); 
    else
        [mcmcInfo, F_array, y_array] = get_state_counts_adj(mcmcInfo); 
    end
%     swap_indices = NaN(n_chains,nStates); % deal with label switching
    
    mcmcInfo = update_v_gibbs(mcmcInfo, F_array, y_array);
    
    % deal with label switches
%     for c = 1:n_chains
%         si = swap_indices(c,:);
%         mcmcInfo.A_curr(:,:,c) = mcmcInfo.A_curr(si,si,c);        
%         mcmcInfo.pi0_curr(c,:) = mcmcInfo.pi0_curr(c,si);
%         
%         sc = mcmcInfo.sample_chains(:,c,:);
%         sc_new = sc;
%         for k = 1:nStates
%             sc_new(sc==k) = si(k);
%         end
%         mcmcInfo.sample_chains(:,c,:) = sc_new; 
%         
%         if update_flag
%             mcmcInfo.pi0_inf_array(update_index,:,c) = mcmcInfo.pi0_inf_array(update_index,si,c);
%             mcmcInfo.v_inf_array(update_index,:,c) = mcmcInfo.v_inf_array(update_index,si,c);
%             mcmcInfo.A_inf_array(:,:,update_index,c) = mcmcInfo.A_inf_array(si,si,update_index,c);
%         end
%     end
    %% %%%%%%%%%%%%% update noise parameter (sigma) %%%%%%%%%%%%%%%%%%%%%%%
    
    % get predicted fluorescence (using new v values)
    mcmcInfo = predict_fluo_full_v3(mcmcInfo);
    
    % Update sigma
    for c = 1:n_chains
                
        % see: https://discdown.org/flexregression/bayesreg.html                              
        a = (numel(mcmcInfo.observed_fluo)/2 + mcmcInfo.a0);    
        F_diff = reshape(permute(mcmcInfo.sample_fluo(us_factor:us_factor:end,c,:),[1 3 2]) - mcmcInfo.observed_fluo,[],1);       
      
        b_prior_piece = mcmcInfo.b0 + (mcmcInfo.v_curr(c,:)-mcmcInfo.v0(c,:))*inv(mcmcInfo.M0)*(mcmcInfo.v_curr(c,:)-mcmcInfo.v0(c,:))';
        b = (F_diff'*F_diff / 2 + b_prior_piece);

        mcmcInfo.sigma_curr(c) = sqrt(1./gamrnd(a,1./b));%mcmcInfo.trueParams.sigma
    
        if update_flag
            mcmcInfo.sigma_inf_array(update_index,c) = mcmcInfo.sigma_curr(c);
        end
    end   
    
    if update_flag && mcmcInfo.save_trace_results
        mcmcInfo.sample_fluo_inf_array(:,:,:,update_index) = mcmcInfo.sample_fluo;
        mcmcInfo.sample_states_inf_array(:,:,:,update_index) = mcmcInfo.sample_chains;
    end