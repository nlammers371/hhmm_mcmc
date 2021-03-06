function mcmcInfo = update_hmm_parameters_par(mcmcInfo)    

    % extrace parameters
    nStates = mcmcInfo.nStates;
    seq_length = mcmcInfo.seq_length;
    n_traces = mcmcInfo.n_traces;
    n_chains = mcmcInfo.n_chains;
    coeff_MS2 = mcmcInfo.coeff_MS2;        
    
    %% %%%%%%%%%% update transition matrix (A) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A_counts = mcmcInfo.transition_count_array;    
    for n = 1:n_chains
        mcmcInfo.A_curr(:,:,n) = sample_A_dirichlet(mcmcInfo.A_alpha, A_counts(:,:,n));
    end       
    
    % update pi0
    for n = 1:mcmcInfo.n_chains
        [V, D] = eig(mcmcInfo.A_curr(:,:,n));
        [~, mi] = max(real(diag(D)));
        mcmcInfo.pi0_curr(n,:) = V(:,mi)/sum(V(:,mi));
    end
    
    update_flag = mod(mcmcInfo.step,mcmcInfo.update_increment) == 0;
    update_index = mcmcInfo.step/mcmcInfo.update_increment + 1;
    if update_flag
        mcmcInfo.A_inf_array(:,:,update_index,:) = mcmcInfo.A_curr;
        mcmcInfo.pi0_inf_array(:,:,update_index,:) = mcmcInfo.pi0_curr;
    end
    
    %% %%%%%%%%%%%%% update emission vector(V) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update V   
    
    % generate F count arrays
    F_array = zeros(seq_length*n_traces,n_chains,nStates);    
    y_array = NaN(seq_length*n_traces,n_chains);    
    
    for n = 1:n_traces
        ind1 = (n-1)*seq_length+1;
        ind2 = n*seq_length;
        % record observed fluo
        y_array(ind1:ind2,:) = repmat(mcmcInfo.observed_fluo(:,n),1,n_chains);
        for m = 1:nStates
            % record counts
            state_counts = convn(coeff_MS2,mcmcInfo.sample_chains(:,:,n)==m,'full');            
            F_array(ind1:ind2,:,m) = state_counts(1:end-length(coeff_MS2)+1,:);                        
        end
    end
       
    % update for each chain
    for n = 1:n_chains
        F_temp = permute(F_array(:,n,:),[1 3 2]);
        M = F_temp'*F_temp + 1e-4;    
        b = F_temp'*y_array(:,n);

        % calculate mean and variance
        v_mean = M\b;
        v_cov_mat = mcmcInfo.sigma_curr(n)^2 * inv(M);

        % sample
        mcmcInfo.v_curr(n,:) = mvnrnd(v_mean, v_cov_mat);             
    end
    if update_flag
        mcmcInfo.v_inf_array(:,:,update_index) = mcmcInfo.v_curr;   
    end
    %% %%%%%%%%%%%%% update noise parameter (sigma) %%%%%%%%%%%%%%%%%%%%%%%

    % get predicted fluorescence
    mcmcInfo = predict_fluo_full_par(mcmcInfo);
    
    for n = 1:n_chains        
        
        % Update sigma
        a = numel(mcmcInfo.observed_fluo)/2;
        F_diff = reshape(permute(mcmcInfo.sample_fluo(:,n,:),[1 3 2]) - mcmcInfo.observed_fluo,[],1);
        b = F_diff'*F_diff / 2;
        
        % fraw sample
        mcmcInfo.sigma_curr(n) = sqrt(1./gamrnd(a,1./b));%mcmcInfo.sigma;%sqrt(mean(F_diff.^2));%
    end
    if update_flag
        mcmcInfo.sigma_inf_array(update_index,:) = mcmcInfo.sigma_curr;
    end
        