function mcmcInfo = update_hmm_parameters_ens(mcmcInfo)    

    % extrace parameters
    nStates = mcmcInfo.nStates;
    seq_length = mcmcInfo.seq_length;
    n_traces = mcmcInfo.n_traces;
    n_chains = mcmcInfo.n_chains;
    coeff_MS2 = mcmcInfo.coeff_MS2;        
    
    %% %%%%%%%%%% update transition matrix (A) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A_counts = mcmcInfo.transition_count_array;        
    mcmcInfo.A_curr = sample_A_dirichlet(mcmcInfo.A_alpha, A_counts);    
    
    % update pi0    
    [V, D] = eig(mcmcInfo.A_curr);
    [~, mi] = max(real(diag(D)));
    mcmcInfo.pi0_curr = V(:,mi)/sum(V(:,mi));    
    
    update_flag = mod(mcmcInfo.step,mcmcInfo.update_increment) == 0;
    update_index = mcmcInfo.step/mcmcInfo.update_increment + 1;
    if update_flag
        mcmcInfo.A_inf_array(:,:,update_index) = mcmcInfo.A_curr;
        mcmcInfo.pi0_inf_array(:,:,update_index) = mcmcInfo.pi0_curr;
    end
    
    %% %%%%%%%%%%%%% update emission vector(V) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update V   
    
    % generate F count arrays
    F_array = zeros(seq_length*n_traces*n_chains,nStates);    
    y_array = NaN(seq_length*n_traces*n_chains,1);    
    
    for n = 1:n_traces        
        ind1 = (n-1)*seq_length*n_chains+1;
        ind2 = n*seq_length*n_chains;
        % record observed fluo
        y_array(ind1:ind2,:) = repmat(mcmcInfo.observed_fluo(:,n),n_chains,1);
        for m = 1:nStates
            % record counts
            state_counts = convn(coeff_MS2,mcmcInfo.sample_chains(:,:,n)==m,'full');            
            F_array(ind1:ind2,m) = reshape(state_counts(1:end-length(coeff_MS2)+1,:),[],1);                        
        end
    end
       
    % update for each chain       
    M = (F_array'*F_array + 1e-4) / n_chains;    
    b = (F_array'*y_array) / n_chains;

    % calculate mean and variance
    v_mean = M\b;
    v_cov_mat = mcmcInfo.sigma_curr^2 * inv(M);

    % sample
    mcmcInfo.v_curr = mvnrnd(v_mean, v_cov_mat)';             
    
    if update_flag
        mcmcInfo.v_inf_array(update_index,:) = mcmcInfo.v_curr;   
    end
    %% %%%%%%%%%%%%% update noise parameter (sigma) %%%%%%%%%%%%%%%%%%%%%%%

    % get predicted fluorescence
    mcmcInfo = predict_fluo_full_ens(mcmcInfo);
       
    % Update sigma
    a = numel(mcmcInfo.observed_fluo)/2;
    F_diff = reshape(permute(mcmcInfo.sample_fluo,[1 3 2]) - mcmcInfo.observed_fluo,[],1);
    b = F_diff'*F_diff / 2 / n_chains;

    % fraw sample
    mcmcInfo.sigma_curr = sqrt(1./gamrnd(a,1./b));%mcmcInfo.sigma;%sqrt(mean(F_diff.^2));%
  
    if update_flag
        mcmcInfo.sigma_inf_array(update_index) = mcmcInfo.sigma_curr;
    end
        