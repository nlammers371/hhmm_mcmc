function mcmcInfo = update_hmm_parameters_v3(mcmcInfo)    

    % extrace parameters
    nStates = mcmcInfo.nStates;
    seq_length = mcmcInfo.seq_length;
    n_traces = mcmcInfo.n_traces;
    n_chains = mcmcInfo.n_chains_eff;
    coeff_MS2 = mcmcInfo.coeff_MS2;        
    
    %% %%%%%%%%%% update transition matrix (A) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    update_flag = mod(mcmcInfo.step,mcmcInfo.update_increment) == 0 || mcmcInfo.step == mcmcInfo.n_mcmc_steps;
    mcmcInfo.update_flag = update_flag;
    update_index = ceil(mcmcInfo.step/mcmcInfo.update_increment) + 1;
    A_counts = mcmcInfo.transition_count_array;  
    ref_chain_ids = repelem(find(mcmcInfo.refChainVec),mcmcInfo.n_temps_per_chain);
    for n = 1:n_chains
        T = 1;%mcmcInfo.tempGradVec(n);        
        A_chain = A_counts(:,:,ref_chain_ids(n)).^(1/T);        
        mcmcInfo.A_curr(:,:,n) = sample_A_dirichlet(mcmcInfo.A_alpha(:,:,n), A_chain);    
        
        % update pi0    
        [V, D] = eig(mcmcInfo.A_curr(:,:,n));
        [~, mi] = max(real(diag(D)));
        mcmcInfo.pi0_curr(n,:) = V(:,mi)/sum(V(:,mi));    
        
        if update_flag
            mcmcInfo.A_inf_array(:,:,update_index,n) = mcmcInfo.A_curr(:,:,n);
            mcmcInfo.pi0_inf_array(update_index,:,n) = mcmcInfo.pi0_curr(n,:);
        end
    end
                
    
    %% %%%%%%%%%%%%% update emission vector(V) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update V   
    
    % generate F count arrays
    F_array = zeros(seq_length*n_traces,nStates,n_chains);        
    y_array = NaN(seq_length*n_traces,n_chains);    
        
    for c = 1:n_chains
        for n = 1:n_traces        
            ind1 = (n-1)*seq_length+1;
            ind2 = n*seq_length;
            % record observed fluo
            y_array(ind1:ind2,c) = mcmcInfo.observed_fluo(:,n);
            for m = 1:nStates
                % record counts
                state_counts = convn(coeff_MS2(:,c),mcmcInfo.sample_chains(:,ref_chain_ids(c),n)==m,'full');            
                F_array(ind1:ind2,m,c) = state_counts(1:end-size(coeff_MS2,1)+1,:);                        
            end
        end
    end  

    for c = 1:n_chains
        T = 1;%mcmcInfo.tempGradVec(n); 
        
        M = ((F_array(:,:,c)'*F_array(:,:,c))) + 1e-4;    
        b = ((F_array(:,:,c)'*y_array(:,c)));
                 
        % calculate mean and variance
        v_mean = M\b;                 
        v_cov_mat = T * inv(mcmcInfo.sigma_curr(c)^-2 * M +  mcmcInfo.sigma_curr(c)^-2 *inv(mcmcInfo.M0));
        
        % sample
        mcmcInfo.v_curr(c,:) = mvnrnd(v_mean, v_cov_mat)'; 
        
        if update_flag
            mcmcInfo.v_inf_array(update_index,:,c) = mcmcInfo.v_curr(c,:);   
        end
    end
    %% %%%%%%%%%%%%% update noise parameter (sigma) %%%%%%%%%%%%%%%%%%%%%%%

    % get predicted fluorescence
    mcmcInfo = predict_fluo_full_v3(mcmcInfo);
    
    % Update sigma
    for c = 1:n_chains
      
        T = 1;%mcmcInfo.tempGradVec(n); 
        % see: https://discdown.org/flexregression/bayesreg.html
        a = (numel(mcmcInfo.observed_fluo)/2 + mcmcInfo.a0).^1/T;    
        
        F_diff = reshape(permute(mcmcInfo.sample_fluo(:,ref_chain_ids(c),:),[1 3 2]) - mcmcInfo.observed_fluo,[],1);       
        b_prior_piece = mcmcInfo.b0 + (mcmcInfo.v_curr(c,:)-mcmcInfo.v0(c,:))*inv(mcmcInfo.M0)*(mcmcInfo.v_curr(c,:)-mcmcInfo.v0(c,:))';
        b = (F_diff'*F_diff / 2 + b_prior_piece)*T;
        
%         beta = 1;
        mcmcInfo.sigma_curr(c) = sqrt(1./gamrnd(a,1./b));%mcmcInfo.trueParams.sigma
        if update_flag
            mcmcInfo.sigma_inf_array(update_index,c) = mcmcInfo.sigma_curr(c);
        end
    end   